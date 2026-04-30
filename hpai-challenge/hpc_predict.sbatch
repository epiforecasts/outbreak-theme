#!/bin/bash
#SBATCH --job-name=hpai_predict
#SBATCH --output=hpai_predict_%j.out
#SBATCH --error=hpai_predict_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

cd ~/hpai-challenge
julia --threads=4 predict.jl 2>&1 | tee predict_output.log
