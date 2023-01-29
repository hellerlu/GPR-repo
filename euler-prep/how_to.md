# Instruction on how to run hybrid jobs on Euler

## Installation of gprMax on Euler
1. Follow instruction on [gprMax Guide](http://docs.gprmax.com/en/latest/include_readme.html):
    1. Download or copy the following files/folder from the gprMax directory to Euler:
        - gprMax
        - build
        - tools
        - user_libs
        - conda_env.yml
        - CREDITS
        - LICENSE
        - README.rst
        - setup.py
    2. Within the Euler directory where you placed the above folders/files, create a new conda environment:
        - make sure you edit the `conda_env.yml` such that mpi4py is under dependencies, rather than pip and uncomment (if applicable)
        - Run: `conda env create -f conda_env.yml`
        - activate conda environment: `conda activate grpMax`
    3. Build and install gprMax:
        - Run: `python setup.py build`
        - Run: `python setup.py install`
    
## Batch file layout
The batch file (shell script .sh) that you need to run a hybrid job according to [Euler documentation](https://scicomp.ethz.ch/wiki/Parallel_computing) looks like: 

> #!/bin/bash
>
> #SBATCH --time=04:00:00  
> #SBATCH --ntasks=8  
> #SBATCH --cpus-per-task=6  
> #SBATCH --ntasks-per-socket=2  
> #SBATCH --mem-per-cpu=2048  
> #SBATCH --nodes=2  
> #SBATCH --job-name=gprmax  
> #SBATCH --output=output_1.txt  
>
> module load openmpi/4.1.4
>
> export OMP_NUM_THREADS=6  
> srun --cpu-bind=cores python -m gprMax files/2D_cylinders_clean.in -n 55 --mpi-no-spawn --geometry-fixed

This runs a hybrid job with 48 cores (N=48), with 8 MPI ranks (M=8) and 6 threads per rank (T=6) and using only 2 MPI ranks per socket (S=2). Also, follow these instruction the get optimal results: 

Let's say you want to run a program on N cores with M MPI ranks and T OpenMP threads per MPI rank where N=M×T. It is strongly advisable that:

1. The number of cores on the node (24 in Euler) is divisible by your chosen T, the number of threads per MPI rank, and
2. You match threads and MPI ranks to the sockets of the node (there are two sockets per node in Euler).

You can also use the API of gprMax, rather than terminal commands and create a runner python file, for example `main_gprMax_run.py`:

> filename = "2D_boxes_clean"
>
> import gprMax  
> gprMax.run(f"files/{filename}.in",n=55,mpi_no_spawn=False, geometry_only=False, geometry_fixed=True)

Make sure to change geometry_fixed=False if you have an antenna model. 

To run this, modify the batch file accordingly: 

> #!/bin/bash
>
> #SBATCH --time=04:00:00  
> #SBATCH --ntasks=8  
> #SBATCH --cpus-per-task=6  
> #SBATCH --ntasks-per-socket=2  
> #SBATCH --mem-per-cpu=2048  
> #SBATCH --nodes=2  
> #SBATCH --job-name=gprmax  
> #SBATCH --output=output_1.txt  
>
> module load openmpi/4.1.4
>
> export OMP_NUM_THREADS=6  
> srun --cpu-bind=cores python main_gprMax_run.py

## Directory structure according to my runner files
Always stay within the directory where you've placed all your gprMax files/folders This is your main working directory. In order to run my runner files, please put a folder `files` in the same directory. Within this `files` folder, put a subfolder called `output_files`. In this folder all outputs will be stored according to its definition in the .in file of the models.  
Put all your models (.in) within the `files` folder. 

