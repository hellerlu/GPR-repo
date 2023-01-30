# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:57:35 2022

Exploring the use of different 

@author: e512481
"""

#%% Generate GPR max files

file_gprmax ="""#title: {title}
#domain: 6 1.8 0.01
#dx_dy_dz: 0.01 0.01 0.01
#time_window: 20e-9

#material: 5 0 1 0 ballast
#material: 7 0 1 0 dry_sand
#material: 15 0 1 0 wet_sand
#material: 25 0.01 1 0 wet_clay
#material: 8 0.01 1 0 concrete
#material: 2 0.01 1 0 dry_wood

#waveform: ricker 1 1.0e9 my_ricker
#hertzian_dipole: z 0.5 0.9 0 my_ricker
#rx: 0.52 0.9 0
#src_steps: 0.02 0 0
#rx_steps: 0.02 0 0
"""

dry_ballast_drysand_wetclay = """
#box: 0 0.6 0 6 0.80 0.01 ballast
#box: 0 0.0 0 6 0.55 0.01 dry_sand
#box: 0 0.0 0 6 0.30 0.01 wet_clay
"""

dry_ballast_wetsand_wetclay = """
#box: 0 0.6 0 6 0.80 0.01 ballast
#box: 0 0.0 0 6 0.55 0.01 wet_sand
#box: 0 0.0 0 6 0.30 0.01 wet_clay
"""

sleepers_steel = """
#box: 1.0 0.6 0 1.145 0.8 0.01 pec
#box: 1.7 0.6 0 1.845 0.8 0.01 pec
#box: 2.4 0.6 0 2.545 0.8 0.01 pec
#box: 3.1 0.6 0 3.245 0.8 0.01 pec
#box: 3.8 0.6 0 3.945 0.8 0.01 pec
#box: 4.5 0.6 0 4.645 0.8 0.01 pec
#box: 5.2 0.6 0 5.345 0.8 0.01 pec
"""
sleepers_concrete = """
#box: 1.0 0.7 0 1.17 0.8 0.01 concrete
#box: 1.7 0.7 0 1.87 0.8 0.01 concrete
#box: 2.4 0.7 0 2.57 0.8 0.01 concrete
#box: 3.1 0.7 0 3.27 0.8 0.01 concrete
#box: 3.8 0.7 0 3.97 0.8 0.01 concrete
#box: 4.5 0.7 0 4.67 0.8 0.01 concrete
#box: 5.2 0.7 0 5.37 0.8 0.01 concrete
""" 
sleepers_wood = """
#box: 1.0 0.64 0 1.16 0.8 0.01 dry_wood
#box: 1.7 0.64 0 1.86 0.8 0.01 dry_wood
#box: 2.4 0.64 0 2.56 0.8 0.01 dry_wood
#box: 3.1 0.64 0 3.26 0.8 0.01 dry_wood
#box: 3.8 0.64 0 3.96 0.8 0.01 dry_wood
#box: 4.5 0.64 0 4.66 0.8 0.01 dry_wood
#box: 5.2 0.64 0 5.36 0.8 0.01 dry_wood
"""

with open("BScan_Railway_track_nosleepers_dryballast25cm_drygravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_drysand_wetclay)
with open("BScan_Railway_track_woodsleepers_dryballast25cm_drygravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_drysand_wetclay+sleepers_wood)
with open("BScan_Railway_track_steelsleepers_dryballast25cm_drygravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_drysand_wetclay+sleepers_steel)
with open("BScan_Railway_track_concretesleepers_dryballast25cm_drygravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_drysand_wetclay+sleepers_concrete)
    
with open("BScan_Railway_track_nosleepers_dryballast25cm_wetgravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_wetsand_wetclay)
with open("BScan_Railway_track_woodsleepers_dryballast25cm_wetgravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_wetsand_wetclay+sleepers_wood)
with open("BScan_Railway_track_steelsleepers_dryballast25cm_wetgravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_wetsand_wetclay+sleepers_steel)
with open("BScan_Railway_track_concretesleepers_dryballast25cm_wetgravel25cm.in",'w') as f:
    f.write(file_gprmax + dry_ballast_wetsand_wetclay+sleepers_concrete)
    
    
    
          
          