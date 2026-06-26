#!/bin/bash
# Parse dependencies passed by snakemake and build a valid Slurm --dependency string.

args="$*"

# If no dependencies provided, print nothing
if [[ -z "$args" ]]; then
  echo -n ""
  exit 0
fi

# Some setups pass Submitted batch job <id> into this field; treat as no deps
if [[ "$args" =~ "Submitted batch job" ]]; then
  echo -n ""
  exit 0
fi

# Extract job IDs
ids=$(grep -Eo '[0-9]{1,10}' <<< "$args" | tr '\n' ',' | sed 's/,$//')

# If nothing extracted, print nothing (avoid --dependency=afterok:)
if [[ -z "$ids" ]]; then
  echo -n ""
  exit 0
fi

echo -n "--dependency=afterok:$ids"

