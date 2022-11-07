from tools.plot_source_wave import *


def cylinder_cmd_block(domain_size_y):
    g = f'''
#python: 
from gprMax.input_cmd_funcs import *

data_file = open("input_files/cirList_1.txt",'r')
for line in data_file:
    cir = line.split()
    cylinder(float(cir[0]), float(cir[1]), 0 , float(cir[0]), float(cir[1]), {domain_size_y}, float(cir[2]), 'ballast')

#end_python:
'''
    return g

def antenna_cmd_block(x,y,z,spatial_res,steps):
    g = f'''
#python:
from user_libs.antennas.GSSI import antenna_like_GSSI_1500
antenna_like_GSSI_1500({x}+{steps}*current_model_run, {y}, {z}, resolution={spatial_res})
#end_python:
    '''
    return g

def command(cmd, *parameters):
    """
    ATTENTION: manipulated so it jumps to the next line automatically
    Helper function. Prints the gprMax #<cmd>: <parameters>. None is ignored
    in the output.

    Args:
        cmd (str): the gprMax cmd string to be printed
        *parameters: one or more strings as arguments, any None values are
            ignored

    Returns:
        s (str): the printed string
    """

    # remove Nones
    filtered = filter(lambda x: x is not None, parameters)
    # convert to str
    filtered_str = map(str, filtered)
    # convert to list
    filtered_list = list(filtered_str)
    try:
        s = '#{}: {} \n'.format(cmd, " ".join(filtered_list))
    except TypeError as e:
        # append info about cmd and parameters to the exception:
        if not e.args:
            e.args = ('', )
        additional_info = "Creating cmd = #{} with parameters {} -> {} failed".format(cmd, parameters, filtered_list)
        e.args = e.args + (additional_info,)
        raise e
    return s

