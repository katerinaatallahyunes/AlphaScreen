#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Katerina Atallah-Yunes

The purpose of this script is to automatically run ChimeraX on AlphaFold output and output contact maps within 
a certain PAE and distance restraints that can be changed by the user
"""
import os
import glob
import sys
from chimerax.core.commands import run

def process_files(pdb_path, pkl_path, protein_prediction, pae, distance):
    print(f"Processing files: PDB - {pdb_path}, PKL - {pkl_path}")
    
    # Verify that the files exist
    if not (os.path.isfile(pdb_path) and os.path.isfile(pkl_path)):
        print(f"Error: One or both of the files do not exist: {pdb_path}, {pkl_path}")
        return
    
    # Load the PDB file
    run(f"open {pdb_path}")
    
    # Load the PAE plot
    run(f"alphafold pae {pkl_path}")
    
    # Custom ChimeraX commands using user input and protein prediction
    run(f"color byattribute {pae}")
    run(f"distance {distance}")
    
    # Example command using protein prediction
    run(f"select {protein_prediction}")
    
    # Command to get the contacts within the PAE and distance restraints
    output_file = f"{protein_prediction}_contacts.txt"
    run(f"alphafold contacts /A to /B distance {distance} maxPae {pae} outputFile {output_file}")

    # Additional commands (customize as needed)
    run("display")
    run("bgcolor white")
    
    # Save the session (modify as needed)
    output_session_path = os.path.join(os.path.dirname(pdb_path), f"{protein_prediction}.cxs")
    run(f"save {output_session_path}")
    
    # Close the current session before processing the next one
    run("close session")

def main(base_directory, pae, distance):
    # Ensure the base directory exists
    if not os.path.exists(base_directory):
        print(f"Error: Base directory '{base_directory}' does not exist.")
        sys.exit(1)

    # Iterate over directories in the base directory
    for root, dirs, files in os.walk(base_directory):
        # Search for pdb and json files
        pdb_files = glob.glob(os.path.join(root, '*.pdb'))
        pkl_files = glob.glob(os.path.join(root, '*.json'))

        if not pdb_files or not pkl_files:
            continue
        
        # Process each pair of pdb and json files
        for pdb_file in pdb_files:
            corresponding_pkl_file = pdb_file.replace('.pdb', '.json')
            if corresponding_pkl_file in pkl_files:
                protein_prediction = os.path.basename(root)  # Use the folder name as protein prediction
                pdb_path = pdb_file
                pkl_path = corresponding_pkl_file
                process_files(pdb_path, pkl_path, protein_prediction, pae, distance)

if __name__ == "__main__":
    # Get user inputs for base directory, PAE, and distance
    if len(sys.argv) != 4:
        print("Usage: python automate_chimerax.py <base_directory> <pae> <distance>")
        sys.exit(1)

    base_directory = sys.argv[1]
    pae = sys.argv[2]
    distance = sys.argv[3]

    # Print inputs for debugging
    print(f"Base Directory: {base_directory}")
    print(f"PAE: {pae}, Distance: {distance}")

    # To run this script in ChimeraX, use the following command:
    # chimerax --nogui --script /path/to/this_script.py <base_directory> <pae> <distance>
    main(base_directory, pae, distance)
