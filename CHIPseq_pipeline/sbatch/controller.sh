#!/bin/bash
#SBATCH --job-name=chip_controller
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/ejvalenzuela/Snakemake/CHIPseq/logs/slurm/controller_CHIPseq_CPT_BW.%j.out
#SBATCH --error=/home/ejvalenzuela/Snakemake/CHIPseq/logs/slurm/controller_CHIPseq_CPT_BW.%j.err

# Controller script used to launch the RNA-seq Snakemake workflow on the SLURM cluster.

module load Miniconda3/24.5.0
eval "$(conda shell.bash hook)"
conda activate snakemake

cd /home/ejvalenzuela/Snakemake/CHIPseq

snakemake \
  --snakefile Snakefile_ChipSeq \
  --profile slurm_config \
  --use-conda --conda-frontend conda \
  --jobscript slurm_config/slurm_jobscript.sh \
  --rerun-incomplete --keep-going \
  --latency-wait 60
