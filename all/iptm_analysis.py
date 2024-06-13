#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Katerina Atallah-Yunes

This script combines two functionalities:
1. Processing AlphaFold prediction folders to extract and manipulate metrics.
2. Using the generated metrics to filter results and generate PAE plots.

Part 1 was adapted from Chop Yan Lee's original code.
The original code source by Chop Yan Lee: https://github.com/KatjaLuckLab/AlphaFold_manuscript/blob/main/scripts/calculate_template_independent_metrics.py
"""

from pymol import cmd
import numpy as np
import matplotlib as plt
import pandas as pd
import json, os, pickle, argparse, sys, csv, shutil, subprocess
from collections import defaultdict

# Part 1: AlphaFold prediction processing

class Prediction_folder:
    """Class that stores prediction folder information"""
    def __init__(self, prediction_folder, num_model=5, project_name=None):
        self.prediction_folder = prediction_folder
        self.num_model = num_model
        self.path_to_prediction_folder = os.path.split(self.prediction_folder)[0]
        self.prediction_name = os.path.split(self.prediction_folder)[1]
        self.rank_to_model = {}
        self.model_confidences = {}
        self.fasta_sequence_dict = {'A':'','B':''}
        self.model_instances = {}
        if project_name is not None:
            self.project_name = project_name
        self.predicted = True

    def parse_ranking_debug_file(self):
        if not os.path.exists(os.path.join(self.prediction_folder, 'ranking_debug.json')):
            self.predicted = False
            return
        else:
            with open(os.path.join(self.prediction_folder, 'ranking_debug.json'), 'r') as f:
                data = json.load(f)
            self.rank_to_model = {f'ranked_{i}': model for i, model in enumerate(data.get("order"))}
            sorted_model_confidence = sorted(data.get("iptm+ptm").values(), reverse=True)
            self.model_confidences = {f'ranked_{i}': float(confidence) for i, confidence in enumerate(sorted_model_confidence)}
        
    def parse_prediction_fasta_file(self):
        fasta_path = f'{self.prediction_folder}.fasta'
        with open(fasta_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip() != '']
        chain_id = 0
        for line in lines:
            if line[0] == '>':
                chain = list(self.fasta_sequence_dict)[chain_id]
                chain_id += 1
                continue
            self.fasta_sequence_dict[chain] += line

    def instantiate_predicted_model(self):
        self.model_instances = {f'ranked_{i}': Predicted_model(f'ranked_{i}') for i in range(self.num_model)}

    def assign_model_info(self):
        for model_id, model_inst in self.model_instances.items():
            model_inst.model_confidence = self.model_confidences.get(model_id)
            model_inst.multimer_model = self.rank_to_model.get(model_id)
            model_inst.path_to_model = self.prediction_folder

    def process_all_models(self):
        self.parse_ranking_debug_file()
        self.parse_prediction_fasta_file()
        if self.predicted:
            self.instantiate_predicted_model()
            self.assign_model_info()
            for model_id, model_inst in self.model_instances.items():
                model_inst.get_model_independent_metrics()
    
    def write_out_calculated_metrics(self, project_name=None):
        metrics_out_path = os.path.join(self.path_to_prediction_folder, 'template_indep_info.tsv')
        metrics_columns_dtype = {'project_name': str, 'prediction_name': str, 'chain_A_length': int, 'chain_B_length': int, 'model_id': str, 'model_confidence': float}
    
        if os.path.exists(metrics_out_path):
            metrics_df = pd.read_csv(metrics_out_path, sep='\t', index_col=0)
            metrics_df.reset_index(drop=True, inplace=True)
        else:
            metrics_df = pd.DataFrame(columns=metrics_columns_dtype.keys())
            metrics_df = metrics_df.astype(dtype=metrics_columns_dtype)
    
        if self.project_name is not None:
            common_info = [self.project_name]
        else:
            common_info = []
        common_info += [self.prediction_name, len(self.fasta_sequence_dict.get('A')), len(self.fasta_sequence_dict.get('B'))]
    
        if not self.predicted:
            row = common_info + ['Prediction failed'] + [None] * 3 + [None] * (len(metrics_columns_dtype) - 6)
            metrics_df.loc[len(metrics_df)] = row
        else:
            for model_id, model_inst in self.model_instances.items():
                row = common_info + [model_id, model_inst.model_confidence]
                metrics_df.loc[len(metrics_df)] = row
    
        filtered_df = metrics_df[metrics_df['model_confidence'] >= 0.5]
    
        filtered_out_path = os.path.join(self.path_to_prediction_folder, 'filtered_template_indep_info.tsv')
        filtered_df.to_csv(filtered_out_path, sep='\t', index=False)
        metrics_df.to_csv(metrics_out_path, sep='\t')
        print(f'Calculated metrics saved in {metrics_out_path}!')
        print(f'Filtered metrics saved in {filtered_out_path}!')

class Predicted_model:
    """Class that stores predicted model"""
    def __init__(self, predicted_model):
        self.predicted_model = predicted_model
        self.path_to_model = None
        self.multimer_model = None
        self.chain_coords = None
        self.chain_plddt = None
        self.pickle_data = None
        self.model_confidence = None

    def check_chain_id(self):
        model_path = os.path.join(self.path_to_model, f'{self.predicted_model}.pdb')
        cmd.load(model_path)
        chains = cmd.get_chains(f'{self.predicted_model}')
        if 'C' in chains:
            cmd.alter(f'{self.predicted_model} and chain B', 'chain="tempA"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain C', 'chain="tempB"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain tempA', 'chain="A"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain tempB', 'chain="B"')
            cmd.sort()
            cmd.save(model_path, f'{self.predicted_model}')
        cmd.reinitialize()

    def read_pickle(self):
        multimer_model_pickle = os.path.join(self.path_to_model, f'result_{self.multimer_model}.pkl')
        with open(multimer_model_pickle, 'rb') as f:
            self.pickle_data = pickle.load(f)

    def parse_atm_record(self, line):
        record = defaultdict()
        record['name'] = line[0:6].strip()
        record['atm_no'] = int(line[6:11])
        record['atm_name'] = line[12:16].strip()
        record['atm_alt'] = line[17]
        record['res_name'] = line[17:20].strip()
        record['chain'] = line[21]
        record['res_no'] = int(line[22:26])
        record['insert'] = line[26].strip()
        record['resid'] = line[22:29]
        record['x'] = float(line[30:38])
        record['y'] = float(line[38:46])
        record['z'] = float(line[46:54])
        record['occ'] = float(line[54:60])
        record['B'] = float(line[60:66])
        return record
    
    def read_pdb(self):
        chain_coords, chain_plddt = {}, {}
        model_path = os.path.join(self.path_to_model, f'{self.predicted_model}.pdb')

        with open(model_path, 'r') as file:
            for line in file:
                if not line.startswith('ATOM'):
                    continue
                record = self.parse_atm_record(line)
                if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                    if record['chain'] in [*chain_coords.keys()]:
                        chain_coords[record['chain']].append([record['x'], record['y'], record['z']])
                        chain_plddt[record['chain']].append(record['B'])
                    else:
                        chain_coords[record['chain']] = [[record['x'], record['y'], record['z']]]
                        chain_plddt[record['chain']] = [record['B']]
        
        for chain in chain_coords:
            chain_coords[chain] = np.array(chain_coords[chain])
            chain_plddt[chain] = np.array(chain_plddt[chain])

        self.chain_coords = chain_coords
        self.chain_plddt = chain_plddt

    def parse_ptm_iptm(self):
        self.ptm = float(self.pickle_data['ptm'])
        self.iptm = float(self.pickle_data['iptm'])

    def get_model_independent_metrics(self):
        self.check_chain_id()
        if 'multimer_v2' in self.multimer_model:
            if os.path.exists(os.path.join(self.path_to_model, f'result_{self.multimer_model}.pkl')):
                self.read_pickle()
                self.parse_ptm_iptm()
            else:
                self.iptm, self.ptm = None, None
        else:
            self.iptm, self.ptm = None, None
        self.read_pdb()

# Part 2: PAE plot generation

def generate_pae_plots(filtered_info_path):
    print(f"Filtered metrics path: {filtered_info_path}")
    filtered_metrics = pd.read_csv(filtered_info_path, sep='\t')
    
    for idx, row in filtered_metrics.iterrows():
        model_id = row['model_id']
        model_confidence = row['model_confidence']
        prediction_name = row['prediction_name']
        prediction_path = os.path.join(filtered_info_path, f'{prediction_name}', model_id)
        
        pae_file_path = os.path.join(prediction_path, 'pae.json')
        
        if not os.path.exists(pae_file_path):
            continue
        
        with open(pae_file_path, 'r') as f:
            pae_data = json.load(f)
        
        pae_matrix = np.array(pae_data['pae'])
        
        plt.imshow(pae_matrix, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.title(f'PAE plot for {prediction_name}, Model {model_id}')
        plt.xlabel('Residue Index')
        plt.ylabel('Residue Index')
        
        plot_path = os.path.join(prediction_path, f'PAE_plot_{model_id}.png')
        plt.savefig(plot_path)
        plt.close()
        print(f'PAE plot saved at {plot_path}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AlphaFold Prediction Processing and PAE Plot Generation")
    parser.add_argument("--prediction_folder", required=True, type=str, help="Path to the prediction folder")
    parser.add_argument("--project_name", required=False, type=str, help="Project name (optional)")
    args = parser.parse_args()

    prediction_folder = args.prediction_folder
    project_name = args.project_name
    
    pred_folder = Prediction_folder(prediction_folder, project_name=project_name)
    pred_folder.process_all_models()
    pred_folder.write_out_calculated_metrics()

    filtered_info_path = os.path.join(pred_folder.path_to_prediction_folder, 'filtered_template_indep_info.tsv')
    generate_pae_plots(filtered_info_path)
