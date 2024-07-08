"""
This script is meant to seperate multiple fastas in a single fasta file to individual fasta files for each protein
"""

import os
import sys

def separate_fasta(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(input_file, 'r') as infile:
        fasta_content = infile.read()

    sequences = fasta_content.split('>')[1:]

    for seq in sequences:
        header, sequence = seq.split('\n', 1)
        sequence = sequence.replace('\n', '')
        identifier = header.split('|')[1]  # Extract the part after the first '|'
        output_file = os.path.join(output_dir, f"{identifier}.fasta")
        
        with open(output_file, 'w') as outfile:
            outfile.write(f">{header}\n{sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python separate_fasta.py <input_file> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    separate_fasta(input_file, output_dir)
