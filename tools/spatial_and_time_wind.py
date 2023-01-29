import numpy as np
from tools.plot_source_wave import *

def get_spatial_and_time_wind(waveform_chosen, amp, freq):
    """ Calculates the spatial resolution and the time window required for the simulation

    Input:
        waveform_chosen:    waveform available in gprMax
        amp:                amplitude of waveform
        freq:               center frequency
    Output:
        Terminal information about resolution and time window
    """
    #timewindow for plot, depends on wavefrom. Ricker usually 4e-9 s
    iter_timewindow = 4e-9 #[s]

    #iteration step, according to http://docs.gprmax.com/en/latest/plotting.html 1.926e-12 works fine:
    iteration_step = 1.926e-12

    #Plots information on chosen waveform
    plot_wave_api(waveform_chosen, amp, freq, iter_timewindow, iteration_step, fft=True)

    #speed of light in vacuum
    c = 2.9979245e8 

    # Highest relative permittivity present in model
    er = 8

    # Maximum frequency present in model, see FFT of waveform, usually 3-4x center frequency
    fmax = 4e9

    # Minimum wavelength, according to https://www.youtube.com/watch?v=DxNITb1Yxyk&t=2706s to mitigate numerical dispertion
    wmin = c / (fmax * np.sqrt(er))

    # Maximum spatial resolution (allowing 10 cells per wavelength)
    dmax = wmin / 10

    #Distance between reciever and bottom of domain and back
    dist_reciever_bottom = 2*1.4 #[m]

    # Pulse width according to wavefrom plotted above
    pulse_width = (2.3-0.5)*1e-9 #[s]

    # Calculated timewindow according to https://www.youtube.com/watch?v=DxNITb1Yxyk&t=2706s: 
    time_wind = dist_reciever_bottom/(c / np.sqrt(6)) + pulse_width
    print('----------------------------------------------------------------')
    print('Minimum wavelength: {:g} m'.format(wmin))
    print('Maximum spatial resolution: {:g} m'.format(dmax))
    print('Minimum timewindow {:g} s'.format(time_wind))
    print('----------------------------------------------------------------')