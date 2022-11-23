#!/bin/bash
#SBATCH -n 40
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --job-name=gprmax
#SBATCH --output=output.txt
#SBATCH --tmp=2000

module load mesa/12.0.6 paraview/5.9.1
export OMP_NUM_THREADS=40; xvfb-run -d pvpython gprMax_pvScript.py
