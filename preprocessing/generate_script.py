import os
import argparse

def generate_slurm_scripts(fasta_folder, output_directory, template_file):
    # Read the template file
    with open(template_file, 'r') as file:
        template = file.read()

    # Iterate through FASTA files in the folder
    for filename in os.listdir(fasta_folder):
        if filename.endswith(".fasta"):
            # Extract the job name from the filename (without the .fasta extension)
            job_name = os.path.splitext(os.path.basename(filename))[0]
            
            # Replace placeholders in the template
            slurm_script = template.format(
                job_name=job_name,
                fasta_path=os.path.join(fasta_folder, filename)
            )
            
            # Define the output script filename
            output_script_filename = os.path.join(output_directory, f"{job_name}_script.sh")
            
            # Write the Slurm script to the output file
            with open(output_script_filename, "w") as output_file:
                output_file.write(slurm_script)
            print(f"Generated Slurm script for {filename}: {output_script_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Slurm scripts for AlphaFold jobs")
    parser.add_argument("fasta_folder", help="Path to the folder containing the FASTA files")
    parser.add_argument("output_directory", help="Directory to save the Slurm scripts")
    parser.add_argument("template_file", help="Path to the template .sh file")
    
    args = parser.parse_args()
    generate_slurm_scripts(args.fasta_folder, args.output_directory, args.template_file)