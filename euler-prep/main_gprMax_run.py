###############################################
# USER INPUT:
#---------------------------------------------#
# Filenames:
filenames = [
    # "2D_boxes_clean",
    # "2D_boxes_clean_antenna"
    "2D_cylinders_clean",
    # "2D_cylinders_clean_antenna"
]
# "3D_cylinders_clean_antenna"
# "3D_boxes_clean"
# "3D_boxes_clean_antenna"
# "3D_cylinders_clean"


# Do you want only the geometry?
geom = False
    

# Is there an antenna involved?
antenna = [
    # False,
    # True,
    False,
    # True
]
###############################################

import gprMax
for i,filename in enumerate(filenames):
    gprMax.run(f"files/{filename}.in",n=55,mpi=15,geometry_only=geom,geometry_fixed=not(antenna[i]))
    