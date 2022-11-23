from tools.classes import *
from tools.input_cmd_funcs import *
from tools.spatial_and_time_wind import get_spatial_and_time_wind
import numpy as np
import os

def create_inputfile():
       
    ################################################################################################
    # USER INPUT:
    #----------------------------------------------------------------------------------------------#

    # Title of B-Scan: 
    title = '2D_cylinders_clean_antenna'

    # Geometry file required? If run as a B-Scan WITH antenna DO NOT include it as it creates a file in each step
    geom_req = False

    # Domain size in x direction
    domain_size_x = 1.5

    # 2D (True) or 3D (False)?
    dim = True

    # Spatial resolution chosen:
    spatial_res = [0.002, 0.002, 0.002]

    # Time window chosen:
    time_window_chosen = 25e-9

    # Antenna or ricker hertzian dipole:
    antenna = True
    
    # Height of source above sleeper:
    height_src_sleeper = 0.7 

    # Source and reciever step per iteration
    # NOTE: NEEDS TO BE A MULTIPLE OF SPATIAL RESOLUTION
    steps = 0.02

    # Geometrical properties of layers:
    gravel_height = 0.15
    asphalt_height = 0.15
    pss_height = 0.25
    ballast_height = 0.25

    # Ballast LOD: ['Box','Cylinder']
    ballast_type = 'Cylinder'

    # Sleeper type: ['concrete','steel','wood']
    sleeper_type = 'concrete'
    # Distance between sleepers (mid-mid)
    dist_sleepers = 0.7 #[m]

    # Distance between domain edge and first sleeper:
    dist_dom_sleeper = 0.3

    # Distance of top ballast to top sleeper:
    dist_lookout = 0.03 # if steel is chosen: MAX 0.01

    # Include rails?
    include_rails = False

    # END USER INPUT
    #################################################################################################

    

    # Filename of geometry file:
    geom_filename = f'{title}_geom'

    # Output directory
    try:
        os.mkdir(f'files/output_files/{title}')
    except FileExistsError:
        print("Output folder already exists!")
    output_dir = f'files/output_files/{title}'
    print(output_dir)

    # Top height of sleepers:
    top_height_sleepers = round(gravel_height+asphalt_height+pss_height+ballast_height+dist_lookout,3)
    
    if dim and not antenna:
        domain_size = [domain_size_x,round(top_height_sleepers+height_src_sleeper+0.2,3),spatial_res[2]]
    elif dim and antenna:
        domain_size = [domain_size_x,round(top_height_sleepers+height_src_sleeper+0.2,3),0.2]   # 3rd dimension according to thickness of antenna
    else:
        domain_size = [domain_size_x,round(top_height_sleepers+height_src_sleeper+0.2,3),1.6]   # fits rails

    if not antenna:
        # Waveform according to http://docs.gprmax.com/en/latest/input.html:
        waveform_chosen = 'ricker'

        #Scaling of the maximum amplitude of the waveform
        amp = 1

        # Center frequency of the waveform in [Hz]:
        freq = 1.0e9

        # Source type
        source_type = 'hertzian_dipole'

        # x-Distance between source and reciever
        dist_src_rx = 0.15
    
        get_spatial_and_time_wind(waveform_chosen,amp,freq)
    
        # Start position of source
        src_start = [0.2, round(top_height_sleepers+height_src_sleeper,3), domain_size[2]/2]  # y=0.9 for 0.4m above ballast
        
        # Calculated number of A-Scans for B-Scan:
        n_steps = int(np.floor((domain_size[0] - 2*src_start[0]-dist_src_rx) / steps))
        print("Number of steps required to run through the whole domain: ", n_steps)
        print('----------------------------------------------------------------')
        print('')
    
    else:     
        # Start position of antenna
        src_start = [0.2, round(top_height_sleepers+height_src_sleeper,3),round(domain_size[2]/2-0.178/2,3)]  #0.178 is GSSi_400 antenna thickness 

        # Calculated number of A-Scans for B-Scan:
        n_steps = int(np.floor((domain_size[0] - 2*src_start[0]) / steps))
        print("Number of steps required to run through the whole domain: ", n_steps)
        print('----------------------------------------------------------------')
        print('')

    
    # Create instances
    mat = materials()
    sl = sleepers(sleeper_type,dist_dom_sleeper,dist_sleepers,top_height_sleepers,domain_size)
    rail = rails(top_height_sleepers,domain_size)

    # Adjust height of aggregates if not already adjusted
    ext_cirList = np.loadtxt('files/cirList_1.txt')
    if any(ext_cirList[:,1] < gravel_height+asphalt_height+pss_height):
        ext_cirList[:,1] = ext_cirList[:,1] + round(gravel_height+asphalt_height+pss_height,2)
    np.savetxt('files/cirList_1.txt',ext_cirList)
    
    # Delete file beforehand, to 'overwrite'
    try:
        os.remove(f'files/{title}.in')
    except:
        pass

    # Create new input file
    f = open(f'input_files/{title}.in','w+')
    f.write('## General comands\n')
    f.write(command('title',title))
    f.write(command('domain',domain_size[0],domain_size[1],domain_size[2]))
    f.write(command('dx_dy_dz',spatial_res[0],spatial_res[1],spatial_res[2]))
    f.write(command('time_window',time_window_chosen))
    f.write(command('output_dir',output_dir))
    f.write('\n')
    if not antenna:
        f.write('## Source:\n')
        f.write(command('waveform',waveform_chosen,amp,freq,f'my_{waveform_chosen}'))
        f.write(command(source_type,'z',src_start[0],src_start[1],src_start[2],f'my_{waveform_chosen}'))
        f.write(command('rx',round(src_start[0]+dist_src_rx,2),src_start[1],src_start[2]))
        f.write(command('src_steps',steps,0,0))
        f.write(command('rx_steps',steps,0,0))
    else:
        f.write('## Antenna:\n')
        f.write(antenna_cmd_block(src_start[0],src_start[1],src_start[2],spatial_res[0],steps))
    
    f.write('\n')
    
    if geom_req and antenna:
        f.write(command('geometry_view',0,0,round(domain_size[2]/2-spatial_res[0],3),domain_size[0],domain_size[1],domain_size[2],spatial_res[0],spatial_res[1],spatial_res[2],geom_filename,'n'))
    else:
        f.write(command('geometry_view',0,0,0,domain_size[0],domain_size[1],domain_size[2],spatial_res[0],spatial_res[1],spatial_res[2],geom_filename,'n'))
    f.write('\n')
    f.write('## Materials:\n')
    f.write(mat.write_materials())
    f.write('\n')
    f.write('## Ballast:\n')
    if ballast_type == 'Box':
        f.write(command('box',0,round(gravel_height+asphalt_height+pss_height,2),0,domain_size[0],round(gravel_height+asphalt_height+pss_height+ballast_height,2),domain_size[2],mat.ballast_mix[0]))
    elif ballast_type == 'Cylinder':
        f.write(cylinder_cmd_block(domain_size[2]))
        f.write(command('box',0,round(gravel_height+asphalt_height+pss_height+ballast_height,2),0,domain_size[0],domain_size[1],domain_size[2],'free_space'))
    f.write('\n')
    f.write('## Subgrade:\n')
    f.write(command('fractal_box',0,round(gravel_height+asphalt_height,3),0,domain_size[0],round(gravel_height+asphalt_height+pss_height,2),domain_size[2],1.5,1,1,1,5,mat.pss[0],f'my_{mat.pss[0]}'))
    f.write(command('add_surface_roughness',0,round(gravel_height+asphalt_height+pss_height,3),0,domain_size[0],round(gravel_height+asphalt_height+pss_height,3),domain_size[2],1.5,1,1,round(gravel_height+asphalt_height+pss_height,2)-0.01,round(gravel_height+asphalt_height+pss_height,2)+0.01,f'my_{mat.pss[0]}'))
    f.write(command('fractal_box',0,gravel_height,0,domain_size[0],round(gravel_height+asphalt_height,3),domain_size[2],1.5,1,1,1,1,mat.asphalt[0],f'my_{mat.asphalt[0]}'))
    f.write(command('add_surface_roughness',0,round(gravel_height+asphalt_height,3),0,domain_size[0],round(gravel_height+asphalt_height,3),domain_size[2],1.5,1,1,round(gravel_height+asphalt_height,3)-0.004,round(gravel_height+asphalt_height,3)+0.004,f'my_{mat.asphalt[0]}'))
    f.write(command('fractal_box',0,0,0,domain_size[0],round(gravel_height,3),domain_size[2],1.5,1,1,1,1,mat.gravel[0],f'my_{mat.gravel[0]}'))
    f.write(command('add_surface_roughness',0,round(gravel_height,3),0,domain_size[0],round(gravel_height,3),domain_size[2],1.5,1,1,round(gravel_height,3)-0.004,round(gravel_height,3)+0.004,f'my_{mat.gravel[0]}'))

    f.write('\n')
    f.write('## Sleepers:\n')
    f.write(sl.write_sleepers())

    # according to rail profil 54 E2 or SBB IV
    if include_rails:
        f.write('\n')
        f.write('## Rails:\n')
        f.write(rail.write_rails())
    
    f.close()

    return title, geom_filename, n_steps