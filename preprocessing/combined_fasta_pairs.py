#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
The purpose of this script is to combine two fastas for an AlphaFold screen. This script differs from the normal combined_fasta, as it combines fastas
from a csv that has different protein pairs, not the same protein against the same proteins.

Input: 
    1. CSV file with all the protein pairs
    2. Path to a folder with all the fasta files for the proteins in the csv file
    3. Path to an output folder where the combined fasta files are saved

Output: 
    1. Folder with all combined fasta files located at the path from input #3
'''

import os
import argparse
import pandas as pd

def is_fasta_file(file_name):
    return file_name.lower().endswith(".fasta")

def read_fasta_file(file_path, encoding='utf-8'):
    sequences = {}
    try:
        with open(file_path, "r", encoding=encoding) as file:
            current_header = ""
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line
                    sequences[current_header] = ""
                else:
                    sequences[current_header] += line
    except Exception as e:
        print(f"Error processing file '{file_path}': {e}")
    return sequences

def write_fasta_file(file_path, sequences, encoding='utf-8'):
    with open(file_path, "w", encoding=encoding) as file:
        for header, sequence in sequences.items():
            file.write(header + "\n")
            file.write(sequence + "\n")

def main(csv_path, fasta_folder_path, output_fasta_path):
    df = pd.read_csv(csv_path)
    if not os.path.exists(output_fasta_path):
        os.makedirs(output_fasta_path)
    
    for index, row in df.iterrows():
        protein1_fasta = os.path.join(fasta_folder_path, row['uid1'] + '.fasta')
        protein2_fasta = os.path.join(fasta_folder_path, row['uid2'] + '.fasta')
        output_file = os.path.join(output_fasta_path, f"{row['uid1']}_{row['uid2']}.fasta")
        
        if os.path.isfile(protein1_fasta) and os.path.isfile(protein2_fasta):
            print(f"Combining '{protein1_fasta}' and '{protein2_fasta}' into '{output_file}'")

            sequences1 = read_fasta_file(protein1_fasta)
            sequences2 = read_fasta_file(protein2_fasta)

            combined_sequences = {**sequences1, **sequences2}

            write_fasta_file(output_file, combined_sequences)
        else:
            if not os.path.isfile(protein1_fasta):
                print(f"'{protein1_fasta}' does not exist in the folder '{fasta_folder_path}'.")
            if not os.path.isfile(protein2_fasta):
                print(f"'{protein2_fasta}' does not exist in the folder '{fasta_folder_path}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine two FASTA files for AlphaFold screening.")
    parser.add_argument("csv_path", help="Path to the CSV file that has all the protein pairs")
    parser.add_argument("fasta_folder_path", help="Path to the folder with all the FASTA files for the proteins to screen against.")
    parser.add_argument("output_fasta_path", help="Path to the folder where the combined FASTA files will be saved.")
    
    args = parser.parse_args()
    
    main(args.csv_path, args.fasta_folder_path, args.output_fasta_path)
