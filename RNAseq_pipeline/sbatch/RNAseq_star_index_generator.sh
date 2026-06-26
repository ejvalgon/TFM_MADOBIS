#!/bin/bash
#SBATCH --job-name=STAR_index_hg38_RNAseq
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ejvalenzuela/Snakemake/RNAseq/logs/star_index_generator/STAR_index_hg38_RNAseq.%j.out
#SBATCH --error=/home/ejvalenzuela/Snakemake/RNAseq/logs/star_index_generator/STAR_index_hg38_RNAseq.%j.err

# ==============================================================================
# Generate STAR genome index for bulk RNA-seq
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Load conda environment
# ------------------------------------------------------------------------------

module load Miniconda3/24.5.0
eval "$(conda shell.bash hook)"
conda activate star_index

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Reference genome FASTA
GENOME_FASTA="/lustre/scratch/profiling/ejvalenzuela/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Gene annotation GTF
GTF_FILE="/lustre/scratch/profiling/ejvalenzuela/genome/human/Homo_sapiens.GRCh38.115.gtf"

# Output STAR index directory
STAR_INDEX="/lustre/scratch/profiling/ejvalenzuela/genome/human/RNAseq_star_index_human"

# Log directory
LOG_DIR="/home/ejvalenzuela/Snakemake/RNAseq/logs/star_index_generator"

mkdir -p "${STAR_INDEX}" "${LOG_DIR}"

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

# For RNA-seq, this should usually be read length - 1.
# Examples:
#   75 bp reads  -> 74
#   100 bp reads -> 99
#   150 bp reads -> 149
SJDB_OVERHANG=149

# ------------------------------------------------------------------------------
# Generate STAR index
# ------------------------------------------------------------------------------

STAR \
  --runMode genomeGenerate \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF_FILE}" \
  --sjdbOverhang "${SJDB_OVERHANG}"

echo "STAR RNA-seq index generated successfully:"
echo "${STAR_INDEX}"
