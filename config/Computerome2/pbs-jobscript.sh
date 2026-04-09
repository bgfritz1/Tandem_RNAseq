#!/bin/sh
# properties = {properties}
conda deactivate 2>/dev/null || true
module load tools snakemake/9.16.3 apptainer/1.4.0 miniconda3/24.9.2
{exec_job}
