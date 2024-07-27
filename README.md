# AlphaScreen

**At the moment, designed for AlphaFold2**

AlphaScreen is a package that allows for streamlining pre- and post-processing analysis of AlphaFold screens. It can be utilized through a single script that does the full analysis. Individual scripts to simply extract metrics, generate PAE plots, and more are also available for individual use. 

**Preprocessing**
Before running any AlphaFold prediction, there are two input files that you need:
  1. FASTA file with the sequences you are predicting
  2. Script to submit the prediction to HPC

All files and scripts to help create input files are located in the `preprocessing` folder. The files in this folder and their descriptions are:

FASTA files / scripts:

  `seperate_fasta.py` - This script seperate multiple fastas in a single fasta file to individual fasta files for each protein.
  
  `chop_fasta.py` - This script chops the full length sequence of a protein into a shorter sequence based on user-input indices.
  
  `combined_fasta.py` - This script combines two fasta files into a single fasta file.    
  
  `combined_fasta_pairs.py` - This script combines two fastas for an AlphaFold screen. This script differs from the normal combined_fasta, as it combines fastas
  from a csv that has different protein pairs, not the same protein against the same proteins.

Script files / scripts:

  `script_template.sh` - This is what a typical AlphaFold script looks like. It **must be modified** based on the user's preferences to work properly. It is one of the inputs required for the `generate_script.sh` script.

  `generate_script.py` - This script uses `script_template.sh` and generates the slurm batch file based, replacing the job name and fasta file based on the fasta file name.
