#!/bin/bash
module load Miniconda3/24.5.0
eval "$(conda shell.bash hook)"
conda activate snakemake

{exec_job}