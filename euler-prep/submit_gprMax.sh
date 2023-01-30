#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-socket=2
#SBATCH --mem-per-cpu=1024
#SBATCH --nodes=6
#SBATCH --job-name=gprmax1
#SBATCH --output=output_1.txt

module load openmpi/4.1.4

export OMP_NUM_THREADS=12
srun --cpu-bind=cores python -m gprMax files/2D_boxes_clean_box_07ghz_wat.in -n 55 --mpi-no-spawn --geometry-fixed

