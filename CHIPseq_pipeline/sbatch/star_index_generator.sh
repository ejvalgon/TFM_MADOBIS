#!/bin/bash
#SBATCH --job-name=STAR_index_spikein
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ejvalenzuela/Snakemake/logs/star_index_generator/STAR_index_human_mouse.%j.out
#SBATCH --error=/home/ejvalenzuela/Snakemake/logs/star_index_generator/STAR_index_human_mouse.%j.err

# ===============================
# Cargar conda y activar el entorno
module load Miniconda3/24.5.0
eval "$(conda shell.bash hook)"
conda activate star_index

# ===============================
# RUTAS

# FASTA combinado humano + raton con prefijos hs_ y mm_
GENOME_FASTA="/lustre/scratch/profiling/ejvalenzuela/genome/human_mouse/human_mouse.fa"

# Directorio donde se guardará el índice STAR combinado
STAR_INDEX="/lustre/scratch/profiling/ejvalenzuela/genome/human_mouse/star_index_human_mouse"

# Crear el directorio del índice si no existe
mkdir -p "${STAR_INDEX}"

# ===============================
# GENERAR ÍNDICE STAR

STAR --runMode genomeGenerate \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${GENOME_FASTA}"
