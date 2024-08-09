#!/bin/bash
#SBATCH --job-name=
#SBATCH --cpus-per-task=
#SBATCH --partition=gpu
#SBATCH --time=
#SBATCH --gres=
#SBATCH --mem=
#SBATCH --mail-user=
#SBATCH --mail-type=

module load 

run_singularity \
--model_preset= \
--fasta_paths= \
--max_template_date= \
--use_precomputed_msas= \
--output_dir= \
--num_multimer_predictions_per_model= \ 
