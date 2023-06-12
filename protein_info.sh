#!/bin/bash

#SBATCH --job-name=Archaea_uniref50_cov50
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G            # --mem=1024G will send to bigmem /// 300G for Uniref50
#SBATCH --output=slurm_output/UR50_comp_output%A.out
#SBATCH --error=slurm_output/UR50_comp_error%A.err


python3 get_protein_id_info_from_subgraphs.py
