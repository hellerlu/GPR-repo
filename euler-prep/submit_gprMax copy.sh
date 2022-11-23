#!/bin/bash
#SBATCH --ntasks=96
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --job-name=gprmax
#SBATCH --output=output_1.txt
#SBATCH --tmp=2000

module load openmpi/4.1.4

export OMP_NUM_THREADS=6
mpirun -n 16 python -m gprMax files/2D_cylinders_clean.in -n 55 --mpi-no-spawn --geometry-fixed
