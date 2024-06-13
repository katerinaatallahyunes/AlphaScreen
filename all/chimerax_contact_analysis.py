#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script filters through PAE / contact data and figures out the most abundant interface within an AlphaFold screen
"""

import os
import re
import csv
from collections import Counter, defaultdict

def extract_values(line):
    """Extracts the values of /A and /B from a line."""
    match = re.search(r'/A:(\d+)\s*/B:(\d+)', line)
    if match:
        a_value = int(match.group(1))
        b_value = int(match.group(2))
        return a_value, b_value
    return None, None

def count_frequencies_in_file(file_path, b_counter, pair_counter, protein_map):
    """Counts frequencies of /A, /B, and (/A, /B) pairs from a single text file."""
    protein_name = os.path.basename(file_path).split('_')[0]  # Extracting protein name from file name
    with open(file_path, 'r') as file:
        for line in file:
            a_value, b_value = extract_values(line)
            if a_value is not None and b_value is not None:
                b_counter[b_value] += 1
                pair_counter[(a_value, b_value)] += 1
                protein_map[b_value].append(protein_name)
            else:
                print(f"Skipped line: {line.strip()}")

def count_frequencies_in_folder(folder_path):
    """Counts frequencies of /A, /B, and (/A, /B) pairs from all text files in a folder."""
    b_counter = Counter()
    pair_counter = Counter()
    protein_map = defaultdict(list)

    if not os.path.isdir(folder_path):
        print("Error: Invalid folder path.")
        return b_counter, pair_counter, protein_map

    files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    if not files:
        print("No files found in the folder.")
        return b_counter, pair_counter, protein_map

    for filename in files:
        file_path = os.path.join(folder_path, filename)
        count_frequencies_in_file(file_path, b_counter, pair_counter, protein_map)

    return b_counter, pair_counter, protein_map

def write_to_csv(output_file, data, protein_map):
    """Writes data to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['B_Value', 'Frequency', 'Proteins'])
        for value, frequency in data.items():
            proteins = ', '.join(protein_map[value])
            writer.writerow([value, frequency, proteins])

def main(folder_path, output_folder):
    b_counter, pair_counter, protein_map = count_frequencies_in_folder(folder_path)

    # Write /B values to CSV
    write_to_csv(os.path.join(output_folder, 'B_values_with_proteins.csv'), b_counter, protein_map)

    # Write (/A, /B) pairs to CSV
    with open(os.path.join(output_folder, 'AB_pairs.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['A_Value', 'B_Value', 'Frequency'])
        for pair, frequency in pair_counter.items():
            writer.writerow([pair[0], pair[1], frequency])

            
if __name__ == "__main__":
    # Replace 'folder_path' with the path to your folder containing text files
    folder_path = '/Users/atallahyuneska/Desktop/test'
    output_folder = '/Users/atallahyuneska/Desktop/test'
    main(folder_path, output_folder)
