#!/bin/bash
#SBATCH --job-name=hpai_validate
#SBATCH --output=hpai_validate_%j.out
#SBATCH --error=hpai_validate_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

cd ~/hpai-challenge
julia --threads=4 validate.jl 2>&1 | tee validate_output.log
