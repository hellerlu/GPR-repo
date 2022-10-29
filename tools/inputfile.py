from tools.classes import *
from tools.input_cmd_funcs import *
import numpy as np

def create_inputfile():
    print('----------------------------------------------------------------')
    # Title of B-Scan: 
    title = '3D_cylinders_clean'

    # Geometry file required? If run as a B-Scan DO NOT include it as it creates a file in each step
    geom_req = True

    # Filename of geometry file:
    geom_filename = f'{title}_geom'

    # Domain size:
    domain_size = [1.5, 1, 1.6]

    # Waveform according to http://docs.gprmax.com/en/latest/input.html:
    waveform_chosen = 'ricker'

    #Scaling of the maximum amplitude of the waveform
    amp = 1

    # Center frequency of the waveform in [Hz]:
    freq = 1.0e9

    # Source type
    source_type = 'hertzian_dipole'

    # Start position of source
    src_start = [0.2, 0.9, domain_size[2]/2]

    # x-Distance between source and reciever
    dist_src_rx = 0.04

    # Source and reciever step per iteration
    # NOTE: NEEDS TO BE A MULTIPLE OF SPATIAL RESOLUTION
    steps = 0.021

    # Calculated number of A-Scans for B-Scan:
    n_steps = int(np.floor((domain_size[0] - 2*src_start[0] - dist_src_rx) / steps))
    print("Number of steps required to run through the whole domain: ", n_steps)
    print('----------------------------------------------------------------')
    print('')

    #timewindow for plot, depends on wavefrom. Ricker usually 4e-9 s
    iter_timewindow = 4e-9 #[s]

    #iteration step, according to http://docs.gprmax.com/en/latest/plotting.html 1.926e-12 works fine:
    iteration_step = 1.926e-12

    #Plots information on chosen waveform
    plot_wave_api(waveform_chosen, amp, freq, iter_timewindow, iteration_step, fft=True)

    #speed of light in vacuum
    c = 2.9979245e8 

    # Highest relative permittivity present in model
    er = 6

    # Maximum frequency present in model, see FFT of waveform, usually 3-4x center frequency
    fmax = 4e9

    # Minimum wavelength, according to https://www.youtube.com/watch?v=DxNITb1Yxyk&t=2706s to mitigate numerical dispertion
    wmin = c / (fmax * np.sqrt(er))

    # Maximum spatial resolution (allowing 10 cells per wavelength)
    dmax = wmin / 10

    #Distance between reciever and bottom of domain and back
    dist_reciever_bottom = 2*0.8 #[m]

    # Pulse width according to wavefrom plotted above
    pulse_width = (2.3-0.5)*1e-9 #[s]

    # Calculated timewindow according to https://www.youtube.com/watch?v=DxNITb1Yxyk&t=2706s: 
    time_wind = dist_reciever_bottom/(c / np.sqrt(6)) + pulse_width
    print('----------------------------------------------------------------')
    print('Minimum wavelength: {:g} m'.format(wmin))
    print('Maximum spatial resolution: {:g} m'.format(dmax))
    print('Minimum timewindow {:g} s'.format(time_wind))
    print('----------------------------------------------------------------')

    # Spatial resolution chosen:
    spatial_res = [0.003, 0.003, 0.003]

    # Time window chosen:
    time_window_chosen = 15e-9

    # Create instances
    mat = materials()

    # Geometrical properties of layers:
    asphalt_height = 0.1
    pss_height = 0.2
    ballast_height = 0.2

    # Ballast LOD: ['Box','Cylinder']
    ballast_type = 'Cylinder'

    # Adjust height of aggregates if not already adjusted
    ext_cirList = np.loadtxt('input_files/cirList_1.txt')
    if any(ext_cirList[:,1] < asphalt_height+pss_height):
        ext_cirList[:,1] = ext_cirList[:,1] + round(asphalt_height+pss_height,2)
    np.savetxt('input_files/cirList_1.txt',ext_cirList)

    ## Sleepers
    # Sleeper type: ['concrete','steel','wood']
    sleeper_type = 'steel'
    # Distance between sleepers (mid-mid)
    dist_sleepers = 0.7 #[m]

    # Distance between domain edge and first sleeper:
    dist_dom_sleeper = 0.3

    # Distance of top ballast to top sleeper:
    dist_lookout = 0.03

    top_height_sleepers = round(asphalt_height+pss_height+ballast_height+dist_lookout,3)
    
    # Create sleeper
    sl = sleepers(sleeper_type,dist_dom_sleeper,dist_sleepers,top_height_sleepers,domain_size)
    
    # Include rails?
    include_rails = True
    rail = rails(top_height_sleepers,domain_size)

    # Delete file beforehand, to 'overwrite'
    try:
        os.remove(f'input_files/{title}.in')
    except:
        pass

    # Create new input file
    f = open(f'input_files/{title}.in','w+')
    f.write('## General comands\n')
    f.write(command('title',title))
    f.write(command('domain',domain_size[0],domain_size[1],domain_size[2]))
    f.write(command('dx_dy_dz',spatial_res[0],spatial_res[1],spatial_res[2]))
    f.write(command('time_window',time_window_chosen))
    f.write(command('waveform',waveform_chosen,amp,freq,f'my_{waveform_chosen}'))
    f.write(command(source_type,'z',src_start[0],src_start[1],src_start[2],f'my_{waveform_chosen}'))
    f.write(command('rx',round(src_start[0]+dist_src_rx,2),src_start[1],src_start[2]))
    f.write(command('src_steps',steps,0,0))
    f.write(command('rx_steps',steps,0,0))
    f.write('\n')
    
    if geom_req:
        f.write(command('geometry_view',0,0,0,domain_size[0],domain_size[1],domain_size[2],spatial_res[0],spatial_res[1],spatial_res[2],geom_filename,'n'))
    f.write('\n')
    f.write('## Materials:\n')
    f.write(mat.write_materials())
    f.write('\n')
    f.write('## Subgrade:\n')
    f.write(command('box',0,0,0,domain_size[0],asphalt_height,domain_size[2],mat.apshalt[0]))
    f.write(command('box',0,asphalt_height,0,domain_size[0],round(asphalt_height+pss_height,2),domain_size[2],mat.pss[0]))

    f.write('\n')
    f.write('## Ballast:\n')
    if ballast_type == 'Box':
        f.write(command('box',0,round(asphalt_height+pss_height,2),0,domain_size[0],round(asphalt_height+pss_height+ballast_height,2),domain_size[2],mat.ballast[0]))
    elif ballast_type == 'Cylinder':
        f.write(cylinder_cmd_block(domain_size[2]))
    f.write(command('box',0,round(asphalt_height+pss_height+ballast_height,2),0,domain_size[0],domain_size[1],domain_size[2],'free_space'))

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