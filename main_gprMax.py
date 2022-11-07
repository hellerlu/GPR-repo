###############################################
# USER INPUT:
#---------------------------------------------#

filename = "3D_cylinders_clean_antenna"
###############################################

import gprMax
gprMax.run(f'input_files/{filename}.in',geometry_only=True)