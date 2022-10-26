# -*- coding: utf-8 -*-
"""
Created on Thu May 13 09:50:37 2021

# Tutorial: https://readgssi.readthedocs.io/en/latest/reading.html
# This package allows to convert gpx (gps data format) to dzg (gssi's propiertary GPS format) https://github.com/iannesbitt/gpx2dzg
@author: e512481
"""
from readgssi import readgssi
from readgssi import plot
import matplotlib.pyplot as plt 
import numpy as np
# import scipy 
# import sys
import GPRfunctions as g
###############################################################################
if __name__ == '__main__' :
    # header, arrs, gps = readgssi.readgssi(r"C:\Users\E512481\Nextcloud\Documents\ETH\ETH_PHD\GPR\Otelfingen_20_12_2021\PROJECT0042012__005.DZT")
    header0, arrs0, gps0 = readgssi.readgssi(r"C:\Users\E512481\Nextcloud\Documents\ETH\ETH_PHD\GPR\Otelfingen_20_12_2021\PROJECT0042012__001.DZT")
    
    header, arrs, gps = readgssi.readgssi(r"C:\Users\E512481\Nextcloud\Documents\ETH\ETH_PHD\GPR\Otelfingen_20_12_2021\PROJECT0042012__003.DZT")
    header1, arrs1, gps1 = readgssi.readgssi(r"C:\Users\E512481\Nextcloud\Documents\ETH\ETH_PHD\GPR\Otelfingen_20_12_2021\PROJECT0042012__006.DZT")
    #%% Dewow filter
    arr = arrs0[0][:,::4]
    fs=np.int32(200/4)
    arrdw = g.dewow(arr, w= fs*50) # 50m dewow
    x = np.arange(0, arr.shape[1]/fs+1/fs ,1/fs)[:-1]
    y = np.arange(0, header0['rhf_range']*(1+1/arr.shape[0]),header0['rhf_range']/arr.shape[0])[:-1]
    
    #%%  Raw plot
    arreg = g.exponential_gain(arrdw, alpha = 0.1, t = np.arange(0, header0['rhf_range']*(1+1/arr.shape[0]),header0['rhf_range']/arr.shape[0])[:-1])#%% plot comparison of raw, dewow and exponential
    fig, ax = plt.subplots(nrows = 1,sharex = True,figsize =(10,6))
    ax = [ax]
    ax[0].pcolormesh(x[100:],y,arr[:,100:], cmap='viridis')
    # ax[1].pcolormesh(x[100:],y,arrdw[:,100:], cmap='viridis')
    # ax[2].pcolormesh(x[100:],y,arreg[:,100:], cmap='viridis')
    
    ax[0].set_ylim(0,20)
    ax2b = ax[0].twinx()
    t0 = 6.972
    v_0 = 0.124
    ax2b.set_ylim((50-t0)*v_0/2,-t0*v_0/2)
    ax2b.set_ylabel('depth (m), v = {}m/ns'.format(v_0))
    
    ax[0].invert_yaxis()
    # ax[1].invert_yaxis()
    # ax[2].invert_yaxis()
    
    
    ax[0].set_title('Raw 25ns Antenna at 0°')
    # ax[1].set_title('Dewow filter, mean 50m')
    # ax[2].set_title('Exponential gain, alpha=0.2')
    
    
    ax[0].set_xlabel('meters')
    # ax[1].set_ylabel('ns')
    ax[0].set_ylabel('ns')
    # ax[2].set_ylabel('ns')
    
    ax[0].set_xlim(14,18.2)
    
    ax[0].set_xlabel('meters')
    fig.savefig("plt/001.DZTraw.png",dpi = 200)
    
    #%% plot comparison of raw, dewow and exponential
    fig, ax = plt.subplots(nrows = 3,sharex = True,figsize =(10,6))
    
    ax[0].pcolormesh(x[100:],y,arr[:,100:], cmap='viridis')
    ax[1].pcolormesh(x[100:],y,arrdw[:,100:], cmap='viridis')
    ax[2].pcolormesh(x[100:],y,arreg[:,100:], cmap='viridis')
    
    ax2b = ax[2].twinx()
    t0 = 6.972
    v_0 = 0.124
    ax2b.set_ylim((50-t0)*v_0/2,-t0*v_0/2)
    ax2b.set_ylabel('depth (m), v = {}m/ns'.format(v_0))
    
    ax[0].invert_yaxis()
    ax[1].invert_yaxis()
    ax[2].invert_yaxis()
    
    
    ax[0].set_title('Raw 25ns Antenna at 0°')
    ax[1].set_title('Dewow filter, mean 50m')
    ax[2].set_title('Exponential gain, alpha=0.2')
    
    
    ax[2].set_xlabel('meters')
    ax[1].set_ylabel('ns')
    ax[0].set_ylabel('ns')
    ax[2].set_ylabel('ns')
    
    fig.savefig("plt/001.DZTraw_dewow_exponetialfilter.png",dpi = 200)
    
    #%% Reflection Coefficients:
    # source of velocities https://gprrental.com/gpr-velocity-table-analysis/
    # meter_per_ns
    v_gpr = {'v_air':0.3,
            'v_metal':10**-12,# 0 --> 100% reflection
            'v_ballast_dry':0.12363652378099924,#estimated
            'v_asphalt_dry':0.17,
            'v_asphalt_wet':0.12,
            'v_concrete_dry':0.13,
            'v_concrete_wet':0.09,
            'v_sand_dry':0.15,
            'v_sand_wet':0.06,
            'v_silt': 0.09,
            'v_clay_dry':0.15,
            'v_clay_wet':0.06,
            'v_organicsoil': 0.04,
            'v_water': 0.033,
            'v_wet_wood': 0.05 # estimated
            } 
        
    R_airsteel = g.reflection_coefficient(v_gpr['v_air'],v_gpr['v_metal'])
    R_airballast = g.reflection_coefficient(v_gpr['v_air'],v_gpr['v_ballast_dry'])
    R_airconcrete = g.reflection_coefficient(v_gpr['v_air'],v_gpr['v_concrete_dry'])
    R_airwetwood = g.reflection_coefficient(v_gpr['v_air'],v_gpr['v_wet_wood'])
    R_ballastsand_dry = g.reflection_coefficient(v_gpr['v_ballast_dry'],v_gpr['v_sand_dry'])
    R_ballastsand_wet = g.reflection_coefficient(v_gpr['v_ballast_dry'],v_gpr['v_sand_wet'])
    # R_ballastclay_dry = reflection_coefficient(v_gpr['ballast'],v_gpr['clay_dry'])
    # R_ballastclay_wet = reflection_coefficient(v_gpr['ballast'],v_gpr['clay_wet'])
    
    #%% Estimated speed from differential reflection at 150m
    # Pythagoras & trigo:
    # singleway_distance = t0*v_air + (t_top-t0)*v_ballast
    # (singleway_distance**2+0.35**2)**0.5 = t0*v_air +(11.8-t0)*v_ballast
    # v_ballast = singleway_distance/(t0*v_air + (t_top-t0))
    # ((t0*v_air)**2-(0.35)**2 + ((11.8-t0)/2*1/(t0*v_air + (t_top-t0))-1)**2*singleway_distance**2 + 2*(t0*v_air)*((11.8-t0)/2*1/(t0*v_air + (t_top-t0)))*singleway_distance) 
    #ax^2+bx+c=0
    
    t0 = 6.972
    v_air = 0.3
    def infer_speedinballast(t0,v_air,t_mid,t_top):
        # t_0 --> time to ground reflection
        # v_air m/s speed in air
        # t_mid --> time to travel between antenna and bottom of ballast on a sleeper
        # t_top --> time to travel between antenna and bottom of ballast between two sleepers
        #   
        # returns speed estimate in ballast
        c = (2*0.35)**2
        b = 2*(t0*v_air)*((t_top-t0)-(t_mid-t0))
        a = ((t_top-t0)**2-(t_mid-t0)**2)
        
        v_ballast = (-b-(b**2-4*a*c)**0.5)/(2*a)
        return v_ballast
    def plot_colormesh_and_speed_inferal(t_mid,t_top,psl):
        
        ps = psl-2.5
        pe = psl+2.5
        
        vv = []
        vs = []
        vs2 = []
        for i in range(200):
            v_ballast= i/600
            singleway_distance = (t0*v_air + (t_top-t0)*v_ballast)/2
            singleway_distance2 = (((t0*v_air +(t_mid-t0)*v_ballast)/2)**2-0.35**2)**0.5
            vv.append(v_ballast)
            vs.append(singleway_distance)
            vs2.append(singleway_distance2)  
        v_ballast = infer_speedinballast(t0,v_air,t_mid,t_top)
        
        fig,ax = plt.subplots( ncols = 2,dpi=200,figsize = (10,6))
        plt.plot(vv,vs, label = 'distance for vertical reflection')
        plt.plot(vv,vs2, label = 'distance for angled reflection 0.35m') 
        plt.scatter(v_ballast,(((t0*v_air +(t_mid-t0)*v_ballast)/2)**2-0.35**2)**0.5,c='b',label= 'effective speed: {}m/ns'.format(round(v_ballast,3)))
        plt.legend()
        plt.xlabel('assumed speed in ballast [m/ns]')
        plt.ylabel('distance traveled [m]')
        
        ax[0].pcolormesh(x[np.argwhere((x>ps) & (x<pe)).flatten()],y,arreg[:,np.argwhere((x>ps) & (x<pe)).flatten()], cmap='viridis')
        ax[0].plot([ps,pe],[t0,t0],label = 'sleeper reflection',c='r')
        ax[0].scatter([psl,psl-0.35],[t_top,t_mid],c='b',label= 'ballast-subbalast reflection')
        
        ax[0].set_ylim(0,20)
        ax[0].set_xlim(ps,pe)
        ax2b = ax[0].twinx()
        ax2b.set_ylim((20-t0)*v_ballast/2,-t0*v_ballast/2)
        ax2b.set_ylabel('depth (m), v = {}m/ns'.format(round(v_ballast,2)))
        ax[0].invert_yaxis()
        ax[0].set_title('Exponential gain, alpha=0.2, position {}m'.format(psl))
        ax[0].set_xlabel('position [m]')
        ax[0].set_ylabel('two way time [ns]')
        ax[0].legend()
        fig.tight_layout()
        return fig,ax
    v_estimated = {}
    position_sleepers = [6.41,21.15,53,79.7,132.28,148,157.95,161.85]
    t_mid_list = [10.2,10.9, 11.2,10.6,11.8,13.6,14.11,15.2]
    t_top_list = [9.15,10.2,10.5,9.9,11.05,12.7,13.2,14.45] 
    for psl, t_mid,t_top in zip(position_sleepers,t_mid_list,t_top_list):
        fig,ax = plot_colormesh_and_speed_inferal(t_mid,t_top,psl)
        v_estimated[psl] = infer_speedinballast(t0,v_air,t_mid,t_top)
        fig.savefig("plt/001.DZT_estimated_wavespeed_ballast_{}.png".format(psl),dpi = 200)
    #%% Plot the estimated speed vs spectrogram
    fig, ax = plt.subplots(nrows =2,sharex = True,figsize =(10,6))
    ax[0].pcolormesh(x[100:],y,arreg[:,100:], cmap='viridis')
    ax[0].plot([x[100],x[-1]],[t0,t0],label = 'sleeper reflection',c='r')
    
    ax[0].set_ylim(0,20)
    ax[0].set_xlim(x[100],x[-1])
    ax2b = ax[0].twinx()
    v_ballast = 0.12
    ax2b.set_ylim((20-t0)*v_ballast/2,-t0*v_ballast/2)
    ax2b.set_ylabel('depth (m), v = {}m/ns'.format(round(v_ballast,2)))
    ax[0].invert_yaxis()
    ax[0].set_title('Exponential gain, alpha=0.2, position {}m'.format(psl))
    ax[1].set_xlabel('position [m]')
    ax[0].set_ylabel('two way time [ns]')
    ax[0].legend()
    ax[1].plot(v_estimated.keys(),v_estimated.values(),label = 'estimated ballast wave speed')
    ax[1].set_ylabel('wave speed ballast (m/ns)')
    fig.tight_layout()
    fig.savefig("plt/001.DZT_estimated_wavespeed_pcolormesh.png",dpi = 200)
    #%% Trying some normalizsation
    plt.figure()
    ar = 0
    plt.plot(arrs[ar][:,:].mean(axis=1))
    array = arrs[ar][:,:].astype("float")
    mean_norm = np.abs((array)).mean(axis=1)
    exp_fit = np.polyfit(np.log(np.arange(511)+1),np.log(mean_norm) ,1)
    expdata = np.exp(exp_fit[1])*np.exp(np.log(np.arange(511)+1)*exp_fit[0])
    plt.plot(np.arange(511), expdata)
             
    #%% Pcolormesh with same normalization as radargram
    import matplotlib.colors as colors
    ar=0
    array = arrs[ar][:,:].astype("float")
    array1 = arrs1[ar][:,::-1].astype("float")
    
    groundreflectionlimit = np.argmax(array.mean(axis=1)) 
    array = array-np.mean(array,axis=1).reshape((-1,1))
    ll = np.min(array[groundreflectionlimit+5:,:])
    ul = np.max(array[groundreflectionlimit+5:,:])
    std = np.std(array[groundreflectionlimit+5:,:])
    groundreflectionlimit  = 0
    gain = 1
    
    
    groundreflectionlimit1 = np.argmax(array1.mean(axis=1)) 
    array1 = array1-np.mean(array1,axis=1).reshape((-1,1))
    # array1 = butter_lowpass_filtfilt(array1, 1, 200)#-np.mean(array1,axis=1).reshape((-1,1)) #200 samples per meter, 1meter cutoff frequency
    
    ll1 = np.min(array[groundreflectionlimit1+5:,:])
    ul1 = np.max(array[groundreflectionlimit1+5:,:])
    std1 = np.std(array[groundreflectionlimit1+5:,:])
    groundreflectionlimit1  = 0
    gain1 = 1
    
    fs = 200/4
    array = array[groundreflectionlimit:,:50000:4]
    array1 = array1[groundreflectionlimit1:,:50000:4]
    
    
    x = np.arange(0, array.shape[1]/fs+1/fs ,1/fs)[:-1]
    y = np.arange(0, header['rhf_range']*(1+1/array.shape[0]),header['rhf_range']/array.shape[0])[:-1]
    
    x1 = np.arange(0, array1.shape[1]/fs+1/fs ,1/fs)[:-1]
    y1 = np.arange(0, header1['rhf_range']*(1+1/array1.shape[0]),header1['rhf_range']/array1.shape[0])[:-1]
    
    fig, ax = plt.subplots(nrows = 2,sharex = True)
    ax[0].pcolormesh(x,y,array, cmap='viridis', clim=(ll, ul),
                         norm=colors.SymLogNorm(linthresh=float(std)/float(gain),
                                                vmin=ll, vmax=ul, base=np.e))
    
    
    a = 5.7 # 
    c = 1
    b = np.sqrt((c/a)**2)
    for i in np.arange(-50,10,1):
        ax[0].plot(np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.33+i*0.66, np.arange(0,25-6.7,0.01)+6.7,'r')
        ax[0].plot(-np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.33+i*0.66, np.arange(0,25-6.7,0.01)+6.7,'r')
    # ax[0].plot(np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.98, np.arange(0,25-6.7,0.01)+6.7,'r')
    # ax[0].plot(-np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.98, np.arange(0,25-6.7,0.01)+6.7,'r')
    ax[0].invert_yaxis()
    
    ax[1].set_xlim(14,22)
    
    ax[1].pcolormesh(x1,y1,array1, cmap='viridis', clim=(ll1, ul1),
                         norm=colors.SymLogNorm(linthresh=float(std1)/float(gain1),
                                                vmin=ll1, vmax=ul1, base=np.e))
    ax[1].invert_yaxis()
    # ax[2].plot(x[3:-2],array1[95,:])#68,:])
    ax[1].set_xlabel('meters')
    ax[1].set_ylabel('ns')
    ax[0].set_ylabel('ns')
    ax[0].set_title('25ns Antenna at 0°')
    ax[1].set_title('25ns Antenna at 90°, lowpass butter filter with cutoff frequency 1[1/m]')
    
    
    
    
    raise Warning('stop')
    #%% plot hyperbola
    plt.figure()
    a = 100/8.1 # maximum distance
    b = 6.7 # distance to first reflection
    plt.plot(np.sqrt((np.arange(0,100,0.05)/a)**2-1*b**2), -np.arange(0,100,0.05)/100*1.4+6.7)
    plt.plot(-np.sqrt((np.arange(0,100,0.05)/a)**2-1*b**2), -np.arange(0,100,0.05)/100*1.4+6.7)
    #%% Concrete Sleeper section
    import matplotlib.colors as colors
    ar=0
    array = arrs[ar][:,:].astype("float")
    array1 = arrs1[ar][:,::-1].astype("float")
    
    groundreflectionlimit = np.argmax(array.mean(axis=1)) 
    array = array -np.mean(array[:,120*200:],axis=1).reshape((-1,1))
    ll = np.min(array[groundreflectionlimit+5:,120*200:])
    ul = np.max(array[groundreflectionlimit+5:,120*200:])
    std = np.std(array[groundreflectionlimit+5:,120*200:])
    groundreflectionlimit  = 0
    gain = 1
    
    
    groundreflectionlimit1 = np.argmax(array1.mean(axis=1)) 
    array1 = array1 -np.mean(array1[:,120*200:],axis=1).reshape((-1,1))
    # array1 = butter_lowpass_filtfilt(array1, 1, 200)#-np.mean(array1,axis=1).reshape((-1,1)) #200 samples per meter, 1meter cutoff frequency
    
    ll1 = np.min(array1[groundreflectionlimit1+5:,120*200:])
    ul1 = np.max(array1[groundreflectionlimit1+5:,120*200:])
    std1 = np.std(array1[groundreflectionlimit1+5:,120*200:])
    groundreflectionlimit1  = 0
    gain1 = 1
    
    fs = 200/4
    array = array[groundreflectionlimit:,:50000:4]
    array1 = array1[groundreflectionlimit1:,:50000:4]
    
    
    x = np.arange(0, array.shape[1]/fs+1/fs ,1/fs)[:-1]
    y = np.arange(0, header['rhf_range']*(1+1/array.shape[0]),header['rhf_range']/array.shape[0])[:-1]
    
    x1 = np.arange(0, array1.shape[1]/fs+1/fs ,1/fs)[:-1]
    y1 = np.arange(0, header1['rhf_range']*(1+1/array1.shape[0]),header1['rhf_range']/array1.shape[0])[:-1]
    
    fig, ax = plt.subplots(nrows = 2,sharex = True)
    ax[0].pcolormesh(x,y,array, cmap='viridis', clim=(ll, ul),
                         norm=colors.SymLogNorm(linthresh=float(std)/float(gain),
                                                vmin=ll, vmax=ul, base=np.e))
    
    
    # a = 5.7 # 
    # c = 1
    # b = np.sqrt((c/a)**2)
    # for i in np.arange(-50,10,1):
    #     ax[0].plot(np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.33+i*0.66, np.arange(0,25-6.7,0.01)+6.7,'r')
    #     ax[0].plot(-np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.33+i*0.66, np.arange(0,25-6.7,0.01)+6.7,'r')
    # # ax[0].plot(np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.98, np.arange(0,25-6.7,0.01)+6.7,'r')
    # # ax[0].plot(-np.sqrt(((np.arange(0,25-6.7,0.01)+c)/a)**2-1*b**2)+17.98, np.arange(0,25-6.7,0.01)+6.7,'r')
    ax[0].invert_yaxis()
    
    ax[1].set_xlim(120,140)
    
    ax[1].pcolormesh(x1,y1,array1, cmap='viridis', clim=(ll1, ul1),
                         norm=colors.SymLogNorm(linthresh=float(std1)/float(gain1),
                                                vmin=ll1, vmax=ul1, base=np.e))
    ax[1].invert_yaxis()
    # ax[2].plot(x[3:-2],array1[95,:])#68,:])
    ax[1].set_xlabel('meters')
    ax[1].set_ylabel('ns')
    ax[0].set_ylabel('ns')
    ax[0].set_title('25ns Antenna at 0°')
    ax[1].set_title('25ns Antenna at 90°, lowpass butter filter with cutoff frequency 1[1/m]')
    
    
    
    #%% Old
    
    plot.radargram(ar =np.divide(array.T,expdata).T,
                   ant=ar, 
                   header=header, 
                   freq=header['antfreq'][ar],
                   verbose=False,
                   figsize=3, 
                   dpi=None, 
                   stack=1, 
                   x='seconds', 
                   z='nanoseconds', 
                   gain=1,
                   colormap= 'viridis',
                   colorbar= True, 
                   noshow=False, 
                   outfile='test_norm', 
                   fmt='png', win=200, 
                   title=True,
                   zero=header['timezero'][ar],
                   zoom=[0,0,0,0], absval=False,
                   showmarks=False)
    
    #%% Plot
    ar = 0
    
    plot.radargram(ar = arrs[ar][:,15000:18000], 
                   ant=ar, 
                   header=header, 
                   freq=header['antfreq'][ar],
                   verbose=False,
                   figsize=3, 
                   dpi=100, 
                   stack=1, 
                   x='meters', 
                   z='nanoseconds', 
                   gain=10,
                   colormap= 'viridis',
                   colorbar= True, 
                   noshow=False, 
                   outfile='test', 
                   fmt='png', win=0, 
                   title=True,
                   zero=header['timezero'][ar],
                   zoom=[0,0,0,0], absval=False, 
                   showmarks=False)
    #%% 21 seconds in one can see the crossing.
    frequency =  header['rhf_sps'] # Hz Sampling frequency in time
    rhb_cdt = header['rhb_cdt']
    gps.iloc[2283] #--> this is the crossing 21 seconds. So the GPS has a 200seconds offset with respect to the radar. 
    # The radar timestamps seem wrong. https://www.google.com/maps/@46.8837293,7.0392841,38m/data=!3m1!1e3
    # -- Ask DB where the 
    
    # For DZT2, concrete bridge at https://www.google.com/maps/place/46%C2%B032'26.5%22N+6%C2%B034'20.0%22E/@46.5425045,6.567852,71m/data=!3m1!1e3!4m13!1m7!3m6!1s0x0:0x0!2zNDbCsDMyJzI2LjUiTiA2wrAzNCcyMC4wIkU!3b1!8m2!3d46.540706!4d6.572232!3m4!1s0x0:0x0!8m2!3d46.540706!4d6.572232
    
    
    #%% Reading the DZT from Edi Meier
    
    # header, arrs, gps = readgssi.readdzt(r"C:\Users\E512481\ownCloud\Documents\ETH\ETH_PHD\GPR\REASSESS\2017-100_651_250M.dt")