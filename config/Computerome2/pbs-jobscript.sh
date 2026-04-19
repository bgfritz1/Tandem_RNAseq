#!/bin/bash
# properties = {properties}

export MAMBA_ROOT_PREFIX=/home/projects/ku_00132/apps/micromamba
module load micromamba/2.5.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

{exec_job}
