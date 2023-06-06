#!/bin/bash

#SBATCH --job-name=Archaea_uniref50_cov50
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G            # --mem=1024G will send to bigmem /// 300G for Uniref50
#SBATCH --output=slurm_output/UR50_comp_output%A.out
#SBATCH --error=slurm_output/UR50_comp_error%A.err

python3 get_connected_components_MA.py ./databases/Uniref/Archaea_uniref50.fasta ./results/mmseq/mmseq_cov/Archaea_Uniref50_covt0.tsv ./results/Archaea_clustered_cov50
echo "now doing cov70"
python3 get_connected_components_MA.py ./databases/Uniref/Archaea_uniref50.fasta ./results/mmseq/mmseq_cov/Archaea_Uniref50_cov70.tsv ./results/Archaea_clustered_cov70
echo "finished"
