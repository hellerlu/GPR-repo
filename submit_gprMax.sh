#!/bin/bash
#SBATCH -n 50
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4096
#SBATCH --job-name=gprmax
#SBATCH --output=output.txt

main_gprMax.py
