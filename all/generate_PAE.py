#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses the iptm_only.py output file to filter through the results, generating PAE plots that fit the user's input threshold
"""

import csv
import os
import subprocess
import shutil

def count_protein_frequencies(file_path, column_name):
    """
    Count the frequencies of protein names in a specified column of a TSV file.

    Args:
    file_path (str): Path to the input TSV file.
    column_name (str): Name of the column containing protein names.

    Returns:
    dict: A dictionary containing protein names and their frequencies in the specified column.
    """
    protein_freq = {}  # Dictionary to store protein frequencies

    with open(file_path, mode='r', newline='') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        
        for row in reader:
            protein_name = row[column_name]
            if protein_name in protein_freq:
                protein_freq[protein_name] += 1
            else:
                protein_freq[protein_name] = 1

    return protein_freq

def find_base_folder_path(current_dir, target_file):
    """
    Find the base folder path containing the target file by traversing up the directory tree.

    Args:
    current_dir (str): The current directory to start the search from.
    target_file (str): The name of the target file to be found.

    Returns:
    str: The absolute path of the base folder containing the target file, or None if not found.
    """
    while True:
        file_path = os.path.join(current_dir, target_file)
        if os.path.isfile(file_path):
            return current_dir
        elif os.path.abspath(current_dir) == os.path.abspath(os.path.join(current_dir, os.pardir)):
            # Reached the root directory
            return None
        else:
            current_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

def copy_fasta_file(base_folder, protein_name):
    """
    Copy the corresponding fasta file for a protein to its respective folder.

    Args:
    base_folder (str): The base folder path containing the protein folders and the fasta file.
    protein_name (str): Name of the protein folder and fasta file.

    Returns:
    bool: True if the fasta file is copied successfully, False otherwise.
    """
    protein_folder_path = os.path.join(base_folder, protein_name)
    fasta_file_src = os.path.join(base_folder, f"{protein_name}.fasta")
    fasta_file_dest = os.path.join(protein_folder_path, f"{protein_name}.fasta")

    if os.path.exists(fasta_file_src) and os.path.isdir(protein_folder_path):
        shutil.copy(fasta_file_src, fasta_file_dest)
        return True
    else:
        return False

def copy_script_file(src_script_path, dest_folder):
    """
    Copy a script file to the destination folder.

    Args:
    src_script_path (str): Path to the source script file.
    dest_folder (str): Path to the destination folder.

    Returns:
    bool: True if the script file is copied successfully, False otherwise.
    """
    if os.path.exists(src_script_path) and os.path.isdir(dest_folder):
        shutil.copy(src_script_path, dest_folder)
        return True
    else:
        return False

def run_script(script_path):
    """
    Run a Python script and capture its output.

    Args:
    script_path (str): Path to the Python script to be executed.

    Returns:
    bool: True if the script runs successfully, False otherwise.
    """
    try:
        result = subprocess.run(['python', script_path], capture_output=True, text=True, check=True)
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running script {script_path}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")
        print(f"Error: {e.stderr}")
        return False

def main():
    file_path = '/Volumes/Untitled/salty_unprocessed/filtered_template_indep_info.tsv'
    column_name = 'prediction_name'  # Name of the column containing protein names
    protein_freq_dict = count_protein_frequencies(file_path, column_name)

    # Set the bound for frequency
    frequency_bound = 3

    # Directory containing the protein folders
    # Find the base folder path dynamically
    current_dir = os.path.dirname(os.path.abspath(file_path))
    base_folder_path = find_base_folder_path(current_dir, 'filtered_template_indep_info.tsv')
    
    if base_folder_path is None:
        print("Error: Could not find the base folder containing the target file.")
        return

    for protein, freq in protein_freq_dict.items():
        if freq >= frequency_bound:
            protein_folder_path = os.path.join(base_folder_path, protein)
            if os.path.exists(protein_folder_path):
                os.chdir(protein_folder_path)  # Change directory to the protein folder
                print(f'Changed directory to: {os.getcwd()}')
                
                # Copy the corresponding fasta file
                if copy_fasta_file(base_folder_path, protein):
                    print(f'Copied fasta file for {protein} to {protein_folder_path}')
                else:
                    print(f'Error copying fasta file for {protein}')

                # Copy the pdockq.py script
                pdockq_script_src = '/path/'
                if copy_script_file(pdockq_script_src, protein_folder_path):
                    print(f'Copied pdockq.py to {protein_folder_path}')
                else:
                    print(f'Error copying pdockq.py to {protein_folder_path}')

                # Copy the plot_AF_all_unrelaxed.py script
                plot_script_src = '/path/'
                if copy_script_file(plot_script_src, protein_folder_path):
                    print(f'Copied plot_AF_all_unrelaxed.py to {protein_folder_path}')
                else:
                    print(f'Error copying plot_AF_all_unrelaxed.py to {protein_folder_path}')

                # Run first Python script
                if not run_script('pdockq.py'):
                    print(f"Error: Script pdockq.py failed for {protein}")
                    continue

                # Run second Python script
                if not run_script('plot_AF_all_unrelaxed.py'):
                    print(f"Error: Script plot_AF_all_unrelaxed.py failed for {protein}")
                    continue

                print(f'Scripts executed for {protein}')
            else:
                print(f'Error: Folder for protein {protein} does not exist!')

if __name__ == "__main__":
    main()

