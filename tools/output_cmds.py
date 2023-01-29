import os
import shutil

def move_output(overwrite,title,geom_filename,n_steps):
    """ Moves .out files into output folder

    Input
        overwrite:      True if existing file should be overwritten
        title:          outputfile title
        geom_filename:  filename of geometry file
        n_steps:        number of files
    """
    # Creating a folder for the files
    try:
        os.mkdir(os.path.join(os.getcwd(),f'output_files/{title}'))
    except FileExistsError:
        if overwrite:
            pass
        else:
            print("Folder already exists, if you want to overwrite its content, set overwrite = True")

    # Delete folder content 
    shutil.rmtree(os.path.join(os.getcwd(),f'output_files/{title}'),ignore_errors=True)

    # Move .out files
    try:
        for i in range(1,n_steps):
            os.rename(f'input_files/{title}{i}.out',f'output_files/{title}/{title}{i}.out')
    except FileNotFoundError:
        print(f"Stopped at: {i}-th .out-file: No files named like this to move.")
    except:
        print("Could not move .out files. Something else went wrong")
    else:
        print("Done with .out files")

    # Move .vti file
    try:
        os.rename(f'input_files/{geom_filename}.vti',f'output_files/{title}/{geom_filename}.vti')
    except FileNotFoundError:
        print(f"No .vit file named like this to move.")
    except:
        print("Could not move .vti file. Something else went wrong")
    else:
        print("Done with .vti file")