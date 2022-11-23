###############################################
# USER INPUT:
#---------------------------------------------#
# Filenames:
filenames = [
    "3D_cylinders_clean_antenna",
    "3D_boxes_clean",
    "3D_boxes_clean_antenna",
    "3D_cylinders_clean"
]

# Do you want only the geometry?
geom = [
    True,
    True,
    True,
    True,
]

# Is there an antenna involved?
antenna = [
    True,
    False,
    True,
    False
]   
###############################################

import gprMax
for i,filename in enumerate(filenames):
    gprMax.run(f'input_files/{filename}.in',geometry_only=geom[i],geometry_fixed=not(antenna[i]))