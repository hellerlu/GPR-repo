# Synthetic GPR data generation for ballast condition assessment



## Introduction
***
This repository contains the code related to my Master thesis at ETH Zürich in HS2022. The main purpose is to generate synthetic GPR data for ballast condition assessment. However, the code can be divded into three sections: 
- Input file creation
    - Varying model approximations
    - Varying simulation parameters
- Running simulation using gprMax
    - Locally 
    - In a HPC environment (ETH Euler Cluster)
- Post-processing analysis
    - Real-world GPR measurements
        - Base-case measurements
        - On-track measurements
    - Simulation output

It features many algorithms to create randomness in input files such as: 
- Random Sequential Adsorption
- Compaction Algortithm

Additionally, many parameters can be adjusted to study their effect on the electromagnetic wave propagation in ballasted track-beds.

The simulation is built upon the [gprMax](https://www.gprmax.com/) software available from Github under: 

https://github.com/gprmax/gprMax


## File structure
***
The GPR-repo has the following file structure:

    gprRepo/
        .vscode/
        data/
            basecase/
            SBB_track/
        euler-prep/
        files/
            output_files/
        received/
        RSA/
            Circles/
            Polygons/
        tools/
        00_main.ipynb
        01_Simulation_analysis.ipynb
        02_RealWorld_Basecase.ipynb
        03_RealWorld_Track.ipynb
        conda_env_gprMax.yml
        conda_env_Repo.yml
        

- `.vscode`: Contains setup and launch json for gprMax dependency (see VS code integration)
- `data`: Contains real-world GPR data from basecase measurments and SBB sample track
- `euler-prep`: Contains samples of batch job code for Euler (read How_to_Euler.md)
- `files`: Folder containing input (.in) and output files (.out) from gprMax simulation
- `received`: Received scripts from supervisors
- `RSA`: Folder containing multiple Python scripts to generate ballast aggregates using RSA and compaction algorithms 
- `tools`: Folder containing several Python scripts that are needed by the 4 Jupyther notebooks
- `00_main`: Main Jupyther notebook to generate input file and running it locally. 
- `01_Simulation_analysis`: Jupyther notebook to compare & analze the simulation output
- `02_RealWorld_Basecase`: Jupyther notebook to compare & analyze basecase measurements
- `03_RealWorld_Track`: Jupyther notebook to compare & analyze on-track measurements
- `conda_env_gprMax`: Conda environment for simulation
- `conda_env_Repo`: Conda environment for analysis


## Setup environment
***
This repo uses Anaconda Navigator to create an evironment. There are two environments:
- gprMax: used to run FDTD simulations and creating input files
- gprRepo: used to run post-processing analysis from simulation and real-world measurement

The loaded packages and dependencies can be seen in the files: 
- *conda_env_gprMax.yml*
- *conda_env_gprRepo.yml*

The can be loaded with the terminal command:
`conda create --name gprRepo --file <file>`

## VS Code integration
***
As mentioned, the simulation is done based on the software [gprMax](https://www.gprmax.com/). The Github repository can be cloned or downloaded. Microsoft Visual Studio code (VS Code) can link two repositories by mapping the other in the PYTHONPATH variable. Two json files need to be created: 
- *launch.json*: Add the following under the configurations:

       {
            "name": "Python-Modul",
            "type": "python",
            "request": "launch",
            "module": "gprMax",
            "justMyCode": true,
            "env":  {
                "PYTHONPATH": "${workspaceFolder}/../gprMax"
            }
        }
- *settings.json": Add the following in the file:

        {  
        "python.analysis.extraPaths": ["${workspaceFolder}/../gprMax"]
        }

Now, all modules from gprMax can be imported as they were in the same repository. 

## How to generate input files
***
Model inputs can be changed and varied in the file `tools/inputfile.py` under the section USER INPUT. If more dielectric properties need to be changed, please do the following:
- Include a variable in the USER INPUT
- Include this variable in the instance creation on line 144
- Go to `tools/classey.py` and include new variable in the material class under the `__init__` method

Everything else is well documented in the `inputfile.py` or in the report.


## Running code on Euler
***
Running the FDTD simulation is integrated for HPC environments. Detailed instructions on how this works on the Euler cluster at ETH can be found in the file: 

`how_to_euler.md` in the *euler-prep* folder

Running Paraview in Server-Client-Mode is difficult and sometimes crashes. But it does work using the instructions in:

`how_to_paraview.md` in the *euler-prep* folder

## General workflow
***
The current general workflow of one synthetic GPR simulation analysis can be descibed as follows:

- **Generate input file** on local machine
    - Use `inputfile.py` and `classes.py` to manipulate model appoximations and material properties
    - Run RSA algorithm to generate new configuration of ballast
    - Follow & run sections in `main.ipynb` to create a new input file
- **Upload .in file to Euler** using [MobaXterm](https://mobaxterm.mobatek.net/)
    - Place .in-File in the folder on Euler
- **Run batch job**
    - According to `How_To_Euler.md`
- **Download all files**
    - Download .out files from cluster using [MobaXterm](https://mobaxterm.mobatek.net/)
    - Download .vti file (geometry)
- **Check geometry**
    - Run paraview locally to inspect models
    - Run paraview in Client-Server Mode according to [EulerDocs](https://scicomp.ethz.ch/wiki/ParaView_Client-Server)
- **Analyse them locally**
    - Use `01_Simulation_analysis`, `02_RealWorld_Basecase`, or `03_RealWorld_Track` to analyze the simulation output with A- or B-Scans. Or do a waveform analysis. 

## Files leftover
***
Most of the .in-files in `files/` have been deleted, especially the ones with different dielectric properties. The ones leftover represent different model approximations and thesis relavant configurations. 

## Contact
***
If any questions arise, please contact Lukas Heller under: 

>lukas.heller(a)windowslive.com


