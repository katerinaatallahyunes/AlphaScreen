#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:18:25 2024

@author: atallahyuneska
"""
import os
import shutil

def copy_fasta_file(base_folder, fasta_folder):
    print(f"Base folder: {base_folder}")
    print(f"FASTA folder: {fasta_folder}")

    # List top-level directories in base_folder
    for item in os.listdir(base_folder):
        item_path = os.path.join(base_folder, item)
        if os.path.isdir(item_path):
            print(f"Analyzing subfolder: {item_path}")
            fasta_file = f"{item}.fasta"
            src_file = os.path.join(fasta_folder, fasta_file)
            dest_file1 = os.path.join(item_path, fasta_file)
            dest_file2 = os.path.join(base_folder, fasta_file)
            print(f"Checking if {src_file} exists")
            if os.path.exists(src_file):
                #print(f"Copying {src_file} to {dest_file}")
                shutil.copyfile(src_file, dest_file1)
                shutil.copyfile(src_file, dest_file2)

                print(f"{fasta_file} has been copied to {item_path}")
            else:
                print(f"{fasta_file} not found in {fasta_folder}")

def main():
    base_folder = 'path/to/folder'
    fasta_folder = 'path/to/folder'
    copy_fasta_file(base_folder, fasta_folder)

if __name__ == "__main__":
    main()
