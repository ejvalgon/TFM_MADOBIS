#!/bin/bash
#SBATCH --job-name=RNAseq_controller
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/ejvalenzuela/Snakemake/RNAseq/logs/slurm/controller_RNAseq_GSE229761.%j.out
#SBATCH --error=/home/ejvalenzuela/Snakemake/RNAseq/logs/slurm/controller_RNAseq_GSE229761.%j.err


# Controller script used to launch the RNA-seq Snakemake workflow on the SLURM cluster.

module load Miniconda3/24.5.0
eval "$(conda shell.bash hook)"
conda activate snakemake

cd /home/ejvalenzuela/Snakemake/RNAseq

snakemake \
  --snakefile Snakefile_RNAseq \
  --profile slurm_config \
  --use-conda --conda-frontend conda \
  --jobscript slurm_config/slurm_jobscript.sh \
  --rerun-incomplete --keep-going \
  --latency-wait 60
