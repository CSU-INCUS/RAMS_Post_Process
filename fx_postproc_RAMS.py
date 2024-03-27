#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 16:37:02 2022

@author: functions used in postproc_RAMS.py
"""

# Import libaries used in the various functions
import h5py
import hdf5plugin
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import re
import datetime
import xarray as xr
import os
from jug import TaskGenerator

def save_rams2D_to_netcdf(h5path,var_names, vars_in, outpath):
 
    # Function that saves extra RAMS variables as compressed NetCDF files
    # Code is based on code from Sean Freeman
    
    # Inputs:
    #   h5file: h5 datafile including full path in str format
    #   vars_in: array of 2D varables to be saved
    #   outpath: pathname where to save netcdf files     
        
    rams_obj = h5py.File(h5path, 'r')

    # xarray needs two dicts: data and coordinates. 
    # Start with 2D coordinates
    coords_out = dict()

    # Pull RAMS model output date from filename
    rams_path_fnonly = h5path.split('/')[-1]
    gridnum = rams_path_fnonly[len(rams_path_fnonly)-4:len(rams_path_fnonly)-3]
    rams_date_str = re.findall("a-[LA]-(.*?)-g"+gridnum+".h5", rams_path_fnonly)[0]
    rams_date = datetime.datetime.strptime(rams_date_str, "%Y-%m-%d-%H%M%S")

    coords_out['time'] = [rams_date,]

 
    coords_out['X'] = np.arange(0, np.shape(vars_in[0])[1], 1)
    coords_out['Y'] = np.arange(0, np.shape(vars_in[0])[0], 1)
    coords_out['lat'] = (['Y', 'X'], rams_obj['GLAT'])
    coords_out['lon'] = (['Y', 'X'], rams_obj['GLON'])

    # Next do the data_out
    data_out = dict()
    # 2D variable, always y by x
    for v in np.arange(0,len(vars_in)):
        data_out[var_names[v]] = (['Y', 'X'], vars_in[v])

    save_data = xr.Dataset(data_vars = data_out, coords=coords_out)

    encoding =  {"zlib": True, "complevel": 6}
    save_filename = 'a-X-'+rams_date_str+'-g'+gridnum+'.nc'
    print('Saving '+save_filename+' to path: '+outpath)
    save_data.to_netcdf(outpath+'/'+save_filename, encoding={var_name: encoding for var_name in list(save_data.variables)}, engine='netcdf4')

    return

def read_var(filename,varname):
    with h5py.File(filename,"r") as f:
        data_out = f[varname][:]
    return data_out

# Need to convert this to a function
def read_head(headfile,h5file):
    
    # Function that reads header files from RAMS
    
    # Inputs:
    #   headfile: header file including full path in str format
    #   h5file: h5 datafile including full path in str format
    
    # Returns:
    #   zmn: height levels for momentum values (i.e., grid box upper and lower levels)
    #   ztn: height levels for thermodynaic values (i.e., grid box centers)
    #   nx:: the number of x points for the domain associated with the h5file
    #   ny: the number of y points for the domain associated with the h5file
    #   npa: the number of surface patches
    
    
    dom_num = h5file[h5file.index('.h5')-1] # Find index of .h5 to determine position showing which nest domain to use

    with open(headfile) as f:
        contents = f.readlines()
        
    idx_zmn = contents.index('__zmn0'+dom_num+'\n')
    nz_m = int(contents[idx_zmn+1])
    zmn = np.zeros(nz_m)
    for i in np.arange(0,nz_m):
        zmn[i] =  float(contents[idx_zmn+2+i])
    
    idx_ztn = contents.index('__ztn0'+dom_num+'\n')
    nz_t = int(contents[idx_ztn+1])
    ztn = np.zeros(nz_t)
    for i in np.arange(0,nz_t):
        ztn[i] =  float(contents[idx_ztn+2+i])
    
    ztop = np.max(ztn) # Model domain top (m)
    
    # Grad the size of the horizontal grid spacing
    idx_dxy = contents.index('__deltaxn\n')
    dxy = float(contents[idx_dxy+1+int(dom_num)].strip())

    idx_npatch = contents.index('__npatch\n')
    npa = int(contents[idx_npatch+2])
    
    idx_ny = contents.index('__nnyp\n')
    idx_nx = contents.index('__nnxp\n')
    ny = np.ones(int(contents[idx_ny+1]))
    nx = np.ones(int(contents[idx_ny+1]))
    for i in np.arange(0,len(ny)):
        nx[i] = int(contents[idx_nx+2+i])
        ny[i] = int(contents[idx_ny+2+i])

    ny_out = ny[int(dom_num)-1]
    nx_out = nx[int(dom_num)-1]

    return zmn, ztn, nx_out, ny_out, dxy, npa 

def convert_soilenergy_to_tempk(so2,w,dryheatcap,npa):

    w = w * 1000 # convert volumentric soil moisture to kg/m^3 by multiplying by density of water (assumed to be 1000)

    # input soil energy (J/m^3)
    sh_l = 4186 #J/(kg*K)
    sh_i = 2093 #J/(kg*K)
    lh_li = 334000 #J/kg

    # Allocate arrays with zero
    frac_liq = np.zeros(np.shape(so2))
    tempk = np.zeros(np.shape(so2))

    # for patch = 0 (water patch, need to do a different calculation)
    if npa == 0:       
        frac_liq[:] = 1
        tempk = (so2 - lh_li) / sh_l + 273.15
    
    else:        
        idx = np.where(so2 <= 0)
        frac_liq[idx] = 0.0
        tempk[idx] = so2[idx] / (sh_i * w[idx] * dryheatcap[idx])  + 273.15
    
        idx = np.where(so2 >= (w * lh_li))
        frac_liq[idx] = 1.0
        tempk[idx] = (so2[idx] - w[idx] * lh_li) / (sh_l * w[idx] + dryheatcap[idx]) + 273.15
    
        idx = np.where((so2 < lh_li) & (so2 > 0.))
        frac_liq[idx] = so2[idx] / (w[idx] * lh_li)
        tempk[idx] = 273.15

    return tempk, frac_liq        

### Calculate surface temperature from RAMS model output fields
def calc_surf_temp(filename):
   
    # Define constants
    extinc_veg = 0.5 # hardcoded, from rams code
    dhc_list = np.array([1465,1407,1344,1273,1214,1177,1319,1227,1177,1151,1088,874])*1000 # Dry eat capacity values for different soil types
   
    # Load model variables
    so_en = read_var(filename,'SOIL_ENERGY')
    so_wa = read_var(filename,'SOIL_WATER')
    veg_temp = read_var(filename,'VEG_TEMP')
    patch_area = read_var(filename,'PATCH_AREA')
    leaf_idx = read_var(filename,'LEAF_CLASS')
    soil_idx = read_var(filename,'SOIL_TEXT')
    veg_frac = read_var(filename,'VEG_FRACAREA') # Only needed for alternate calculate of surface temperature
    
    i_sl = np.shape(so_en)[1]-1 # soil level closest to the surface (-1 since python uses 0-based index)
 
    # Calculate soil temperature based on soil energy 
    tempk = OrderedDict()
    frac_liq = OrderedDict()   
    # Loop through number of patches
    for npa in np.arange(0,np.shape(patch_area)[0]):
        dhc = np.zeros(np.shape(so_en[npa,i_sl,:,:])) # dry heat capacity of soil for different soil types
        for i in np.arange(0,np.shape(dhc)[0]):
            for j in np.arange(0,np.shape(dhc)[1]):
                dhc[i,j] = dhc_list[int(soil_idx[npa,i_sl,i,j]-1)] # Dry Heat Capactity, based on soil type.   
        tempk[npa], frac_liq[npa] = convert_soilenergy_to_tempk(so_en[npa,i_sl,:,:],so_wa[npa,i_sl,:,:],dhc,npa)
    
    # Temperature (K) full, as a linear combination of temperature of individual patches
    tempk_f = np.zeros(np.shape(tempk[0])) 
    for npa in np.arange(0,np.shape(patch_area)[0]):
        tempk_f = tempk_f + patch_area[npa,:,:] * tempk[npa]
               
    # Vegatation Fraction as a function of LEAF_CLASS 
    # This does not consider radiative effects through the canopy for LW, which VEG_FRAC variable does    
    vf_list = [0,0,0,0,0.8,0.8,0.8,0.9,0.75,0.8,0.2,0.6,0.7,0.7,0.8,0.85,0.8,0.8,0.8,0.74,0.9] # from 
    vf = np.zeros(np.shape(so_en[npa,i_sl,:,:]))
    # Loop through x,y, to determine vegation fraction based on leaf class at each model grid point
    for i in np.arange(0,np.shape(vf)[0]):
        for j in np.arange(0,np.shape(vf)[1]):
            vf[i,j] = vf_list[int(leaf_idx[1,i,j])]
    
    # Calculate surface temperature as a combination of soil temperature and vegetation temperature based on the vegetation fraction
    # and the soil fraction (i.e., 1 - vegetation fraction)
    tempk_f3 = np.zeros(np.shape(tempk[0]))
    for npa in np.arange(0,np.shape(patch_area)[0]): # Loop through patches
        tempk_f3 = tempk_f3 + (vf * veg_temp[npa,:,:] + (1-vf) * tempk[npa]) * patch_area[npa,:,:]
        
    # Alternatve calculate of surface temperature using veg_frac variable, which considers the amount of LW radation that makes it through the canopy
    #tempk_f3 = np.zeros(np.shape(tempk[0]))
    #for npa in np.arange(0,3):
    #    tempk_f3 = tempk_f3 + (veg_frac[npa,:,:] * veg_temp[npa,:,:] + (1-veg_frac[npa,:,:]) * tempk[npa]) * patch_area[npa,:,:]

    return tempk_f3
        
### Calculate 10m wind speed from RAMS model output fields
def calc_10m_wind(filename,ztn):

    # Define constants
    vonk = 0.40 # von Karman constant
    cp = 1004 # J/kg/K
    g = 9.8 # m/s

    # Use RAMS model heights (grid midpoints, ztn) to interpolate values to 10 m above ground
    zold = ztn[1] # RAMS height just above the surface
    znew = 10. # 10 m height above the surface
    ztop = np.max(ztn) # Model domain top (m)
   
    # Load model variables
    pi = read_var(filename,'PI')
    topo = read_var(filename,'TOPT')
    prough = read_var(filename,'PATCH_ROUGH')
    patch_frac = read_var(filename,'PATCH_AREA')
    cantemp = read_var(filename,'CAN_TEMP')
    u = read_var(filename,'UP') # u component of wind (m/s)
    v = read_var(filename,'VP') # v component of wind (m/s)
    theta = read_var(filename,'THETA') # Potential temperature (K)
    ustar = read_var(filename,'USTAR') # Read in friction velocity (m/s)
 
    # Calculate wind speed from u and v components
    speed = np.sqrt(u*u+v*v)
     
    speed_10m = np.zeros(np.shape(topo))  
    
    # Save intermediate variables for testing
    #richno_s = np.zeros([int(nx[gi]),int(ny[gi])])
    #factor_s = np.zeros([3,int(nx[gi]),int(ny[gi])])
    #factorpos_s = np.zeros([3,int(nx[gi]),int(ny[gi])])
    #factorneg_s = np.zeros([3,int(nx[gi]),int(ny[gi])])
    #z0_s = np.zeros([3,int(nx[gi]),int(ny[gi])])
    
    # Loop through x,y domain and calculate 10 m wind at each grid point
    for i in np.arange(0,np.shape(speed_10m)[0]):
        for j in np.arange(0,np.shape(speed_10m)[1]):
            
            rtgt = 1.0-(topo[i,j]/ztop) # Topography adjustment        
            sfcpi = 0.5 * (pi[0,i,j]+pi[1,i,j])# Calculate surface RAMS exner function
            zagl = zold*rtgt # Adjust RAMS heights for topography
            
            rwindp = 0 # Initial 10 m wind value

            # Loop through patches
            for p in np.arange(0,np.shape(patch_frac)[0]):
                
                z0 = prough[p,i,j] # Patch roughness length for each grid point
                if p == 0:
                    z0 = 0.001 # For water patch, force roughness length to 0.001 m (found in RAMS code)
                
                spd = np.max([speed[1,i,j],0.25]) # Make minimum speed 0.25 m/s for the second model level / Note second model level is the same as first model level, which is below ground
                cantheta = cantemp[p,i,j]*cp/sfcpi # Convert canopy temperature to canopy potential temperature
                
                # Calculate a variant of the Richardson Number
                richno = g * zagl * (theta[1,i,j] - cantheta) / (theta[1,i,j] * np.power(spd,2.0))
                
                a2 = np.power( vonk / np.log(znew / z0) , 2.0 ) # Ratio of ustar / u (squared), see  https://en.wikipedia.org/wiki/Von_K%C3%A1rm%C3%A1n_constant
    
                rwind_init = np.power(ustar[p,i,j],2.0) / a2 # by dividing by a2, ustar's cancel and left with u at 10m
                
                # Apply factors to the derived wind, which vary if the richardson number is positive or negative
                if richno > 0: # if Richardson nmber is positive (Atmosphere air is warmer than canopy air)
                    factor = (1.0+10.0*richno/np.sqrt(1+5*richno))
                    rwind = np.sqrt(rwind_init * factor)
                    # Save intermediate variables for testing    
                    #factorpos_s[p,i,j] = np.sqrt(factor)
                    #factorneg_s[p,i,j] = np.nan
                else: # if Richardson nmber is negative or 0                    
                    factor = (1.0-10.0*richno/(1.0+75.0*a2*np.sqrt(-znew*richno/z0)))
                    rwind = np.sqrt(rwind_init / factor)
                    # Save intermediate variables for testing    
                    #factorneg_s[p,i,j] = np.sqrt(factor)
                    #factorpos_s[p,i,j] = np.nan

                # Save intermediate variables for testing    
                #factor_s[p,i,j] = np.sqrt(factor)
                #z0_s[p,i,j] = z0
    
                rwind = np.max([ np.min([rwind,speed[1,i,j]]) , 0.0]) # wind is at most the initial value of the surface wind, and must be greater than 0
                rwindp = rwindp + rwind * patch_frac[p,i,j] # sum over patches
    
            speed_10m[i,j] = rwindp

    return speed_10m

@TaskGenerator
def plot_w_itc(h5file,savepath,zt,run_ID):

    cur_time = os.path.split(h5file)[1][9:21]

    # Constants for calculating total integrated condensate
    cp = 1004; # J/kg/K
    rd = 287; # J/kg/K
    p00 = 100000; # Reference Pressure
    alt_des = 5000;

    aid = np.where(np.abs(np.array(zt)-alt_des)==np.min(np.abs(np.array(zt)-alt_des)))[0]

    lat = read_var(h5file,'GLAT')
    lon = read_var(h5file,'GLON')
 
    wp = read_var(h5file,'WP') # u component of wind (m/s)
    nx = np.shape(wp)[2]
    ny = np.shape(wp)[1]
    # Sum hydrometeor mass mixing ratios to calculate total mass mixing ratio
    #rtp = read_var(h5file,'RCP')+read_var(h5file,'RRP')+read_var(h5file,'RPP')+read_var(h5file,'RSP')+read_var(h5file,'RAP')+read_var(h5file,'RHP')+read_var(h5file,'RGP')
    rtp = read_var(h5file,'RTP')-read_var(h5file,'RV')    

    # Load variables needed to calculate density
    th = read_var(h5file,'THETA')
    pi = read_var(h5file,'PI')
    rv = read_var(h5file,'RV')
    
    # Convert RAMS native variables to temperature and pressure
    pres = np.power((pi/cp),cp/rd)*p00
    temp = th*(pi/cp)
    del(th,pi)
    
    # Calculate atmospheric density
    dens = pres/(rd*temp*(1+0.61*rv))
    del(pres,temp,rv)
    
    # Difference in heights (dz)    
    diff_zt_3D = np.tile(np.diff(zt),(int(ny),int(nx),1))
    diff_zt_3D = np.moveaxis(diff_zt_3D,2,0)

    # Calculate integrated condensate
    itc = np.nansum(rtp[1:,:,:]*dens[1:,:,:]*diff_zt_3D,axis=0) # integrated total condensate in kg
    itc_mm = itc/997*1000 # integrated total condensate in mm
    itc_mm[itc_mm<=0] = 0.001
    w_plt = np.nanmax(wp[int(aid):,:,:],axis=0)
    
    # contour and colobar ticks and levels     
    w_lvls1 = np.array([2])
    w_lvls2 = np.array([10])
    w_lvls3 = np.array([30])
    itc_lvls = np.arange(0.01,10.01,0.01) # Adjusted these levels, such that figure shows regions with at least 1 grid box with 0.1 g/kg of condensate
    itc_cbar_ticks = np.log10(np.array([1,5,10]))
    itc_cbar_ticklbls = np.array([1,5,10])
    
    # Make new colorbar to blue (no condensate) to white (condensate)
    from matplotlib.colors import LinearSegmentedColormap
    colorlist=["darkblue", "lightsteelblue", "white"]
    newcmp = LinearSegmentedColormap.from_list('testCmap', colors=colorlist, N=256)
    
    # Scale size of figure based on dimensions of domain
    max_dim = np.max([nx,ny])
    fs_scale = 9
    lw = 1.0
   
    # Plot Figure
    fig,ax = plt.subplots(1,1,figsize=[nx/max_dim*fs_scale,ny/max_dim*fs_scale/1.2]) # Divide y by 1.2 to account for colorbar
    #a = ax.contourf(np.log10(itc_mm),levels=np.log10(itc_lvls),extend='both',cmap=plt.cm.Blues)
    ax.grid()
    #a = ax.contourf(np.log10(itc_mm),levels=np.log10(itc_lvls),cmap=newcmp,extend='both')
    #b = ax.contour(w_plt,levels=w_lvls1,colors='gold',linewidths=lw)
    #b = ax.contour(w_plt,levels=w_lvls2,colors='m',linewidths=lw)
    #b = ax.contour(w_plt,levels=w_lvls3,colors='limegreen',linewidths=lw)
    a = ax.contourf(lon,lat,np.log10(itc_mm),levels=np.log10(itc_lvls),cmap=newcmp,extend='both')
    b = ax.contour(lon,lat,w_plt,levels=w_lvls1,colors='gold',linewidths=lw)
    b = ax.contour(lon,lat,w_plt,levels=w_lvls2,colors='m',linewidths=lw)
    b = ax.contour(lon,lat,w_plt,levels=w_lvls3,colors='limegreen',linewidths=lw)
    cbar = plt.colorbar(a,ax=ax,ticks=itc_cbar_ticks)
    cbar.ax.set_yticklabels(itc_cbar_ticklbls)
    cbar.ax.set_ylabel('Integrated Total Condensate (mm)')
    ax.set_title('Integrated Total Condensate (Shaded) @ '+cur_time+' \n Maximum W above 5 km (Contoured) \n (Gold, purple and green contours: 2, 10, and 30 m s$^{-1}$ , respectively)')
    ax.set_xlabel(run_ID)
    plt.tight_layout()
    fig.savefig(savepath+'/'+'WMax_IntCond_'+cur_time+'.png')
    plt.close(fig)
       
    return
