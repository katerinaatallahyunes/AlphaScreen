#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task=
#SBATCH --partition=gpu
#SBATCH --time=
#SBATCH --gres=
#SBATCH --mem=
#SBATCH --mail-user=
#SBATCH --mail-type=

module load 

run \
--model_preset= \
--fasta_paths={fasta_path} \
--max_template_date= \
--output_dir= \
--num_multimer_predictions_per_model= 
