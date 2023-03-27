"""
Author: Peter J. Marinescu
Script to plot plan views from RAMS model output for 2D or 3D variables

"""

# Import python libraries
import numpy as np
import matplotlib.pyplot as plt
import glob
import rams_tools
import os
import h5py
import hdf5plugin
from collections import OrderedDict
import math
from jug import TaskGenerator

# Path to RAMS .h5 and head.txt files
#path = '/nobackup/ldgrant/pmarines/PROD/ARG1.2-R/G3/out_30s/'
path = '/nobackup/pmarines/PROD/DRC1.1-R/G3-C/out_30s/'

# Path to save plots
#savepath = '/nobackup/ldgrant/pmarines/PROD/ARG1.2-R/G3/out_30s/Plots/'
savepath = '/nobackup/pmarines/PROD/DRC1.1-R/G3-C/out_30s/Plots/'

# Specify array of variables to plots
# ****** Make sure to check and make sure variable is present in clvls, cmaps, and scale dictionaries below; If not there add it
vars = ['VP','UP','WP'] # RAMS model variables
#vars = ['RTP','RLONTOP']

# Specify altitudes of levels to plots
lvls = [50,5000,10000,15000] # m

# Specify grid number to plot
gridnum = '3'
filetype = 'L' # Analysis (A) or Lite (L) files

# Specify grid point bounds for zooming in
# leave x0 = -1 to not use this feature
x0 = -1; x1 = 0;
y0 = 0; y1 = 0;
sname = ''; # ability to add an additional phrase to filename associated with zooming in

# Use Lat/lons or grid points for x,y axes
use_latlon = 1; #1 for lat lons, 0 for grid points

# Specify either an array or single integer to specify contour levels or number of levels for each variable
clvls = OrderedDict()
clvls['WP'] = np.arange(-40,40.1,0.5)
clvls['UP'] = np.arange(-40,40.1,0.5)
clvls['VP'] = np.arange(-40,40.1,0.5)
clvls['PCPVR'] = np.arange(0,100.1,1.0)
clvls['RTP'] = np.arange(0,20,0.5)
clvls['RV'] = np.arange(0,20,0.5)
clvls['RLONTOP'] = np.arange(60,300.1,5)

cmaps = OrderedDict()
cmaps['WP'] = plt.cm.bwr
cmaps['UP'] = plt.cm.bwr
cmaps['VP'] = plt.cm.bwr
cmaps['PCPVR'] = plt.cm.viridis
cmaps['RTP'] = plt.cm.viridis
cmaps['RV'] = plt.cm.viridis
cmaps['RLONTOP'] = plt.cm.viridis

scale = OrderedDict()
scale['WP'] = 1 #m/s
scale['UP'] = 1 #m/s
scale['VP'] = 1 #m/s
scale['PCPVR'] = 1 / 997 * 1000 * 3600 #mm/hr
scale['RTP'] = 1000 #g/kg
scale['RV'] = 1000 #g/kg
scale['RLONTOP'] = 1 #W/m2

# Specify details of a potential new domain
plt_mult_dom = 0; # 1 for plotting multiple available analysis domains, 0 for plotting just specified domain above
plt_new_dom = 0; # 0 for no, 1 for yes; only works when plt_mult_dom is on
lat_p = OrderedDict() ; lon_p = OrderedDict()
name = 'new'; sname3 = '' # Don't Change these parameters
cen_lat = -31.75; cen_lon = -64.5
nx = 2630; ny = 2780;
dx = 100 ; dy = 100 #m

cen_lat = -4.0; cen_lon = 23.7
nx = 1424; ny = 1208;
dx = 400 ; dy = 400 #m

cen_lat = -4.3; cen_lon = 24.9
nx = 1880; ny = 2670;
dx = 100 ; dy = 100 #m

#### Define Function for plotting single domain (for parallization)
@TaskGenerator
def plot_single_dom(var, znow, hefile, h5file, x0,x1,y0,y1, scale, clvls, cmaps, tlvls, savepath, sname, gridnum):

    # Get RAMS vertical coodinate levels
    zcoords = np.array(rams_tools.calc_zcoords(hefile))
    dzcoords = np.diff(zcoords)
    # Get the vertical coordinate index associated with the specified altitudes
    zid = np.where(np.abs(zcoords-znow) == np.min(np.abs(zcoords-znow)))[0][0]
    print(zid)

    # convert altitude to string for filenames
    lvl = int(np.array(zcoords[zid]))
    if lvl < 10:
       lvlstr = '0000'+str(lvl)
    elif lvl < 100:
       lvlstr = '000'+str(lvl)
    elif lvl < 1000:
       lvlstr = '00'+str(lvl)
    elif lvl < 10000:
       lvlstr = '0'+str(lvl)
    else:
       lvlstr = str(lvl)

    # Grab current time from rams filename
    cur_time = os.path.split(h5file)[1][9:21]
    print(cur_time)

    # Open RAMS file
    rams_file = h5py.File(h5file, 'r')

    # print variables to test code
    print(np.shape(rams_file[var]))
    print(len(np.shape(rams_file[var])))

    # Use zoomed in region; else use entire domain
    if x0 >= 0:
        lats = np.array(rams_file['GLAT'][x0:x1,y0:y1])
        lons = np.array(rams_file['GLON'][x0:x1,y0:y1])
        topt = np.array(rams_file['TOPT'][x0:x1,y0:y1])
        if len(np.shape(rams_file[var])) == 2:
            pltvar = np.array(rams_file[var][x0:x1,y0:y1])
        else:
            pltvar = np.array(rams_file[var][zid,x0:x1,y0:y1])
    else:
        lats = np.array(rams_file['GLAT'][:,:])
        lons = np.array(rams_file['GLON'][:,:])
        topt = np.array(rams_file['TOPT'][:,:])
        if len(np.shape(rams_file[var])) == 2:
            pltvar = np.array(rams_file[var][:,:])
        else:
            pltvar = np.array(rams_file[var][zid,:,:])

    # use grid points for x,y axes
    if use_latlon == 0:
        lons = np.arange(0,np.shape(topt)[1])
        lats = np.arange(0,np.shape(topt)[0])
        sname2 = 'GridPt'
    else:
        sname2 = 'LatLon'

    fsx = 8;
    fsy = fsx*np.shape(topt)[0]/np.shape(topt)[1]/1.2
    # Make Figure
    fig,ax = plt.subplots(1,1,figsize=(fsx,fsy))

    a = ax.contourf(lons,lats,pltvar*scale[var],levels=clvls[var],cmap=cmaps[var],extend='both')
    b = ax.contour(lons,lats,topt,levels=tlvls,colors='k')

    plt.colorbar(a,ax=ax)
    ax.set_title(var+' @ '+str(np.round(zcoords[zid]/1000,1))+' km:'+cur_time)
    ax.grid()
    plt.tight_layout()

    plt.savefig(savepath+sname2+sname+'G'+gridnum+'_'+var+'_'+lvlstr+'_'+cur_time+'.png')
    plt.close(fig)

    return

#### Define Functions for creating and plotting new grid ideas

def domain_box_from_center(cen_lat,cen_lon,dy,dx):

    factors = [[1,-1],[1,1],[-1,1],[-1,-1]]
    Nlon = np.zeros(4)
    Nlat = np.zeros(4)
    # loop through four corners
    for i in np.arange(0,4):
        dy_adj = dy * factors[i][0]
        dx_adj = dx * factors[i][1]
#        print(factors[i][0],factors[i][1])
#        print(dx,dy)
#        print(dx_adj,dy_adj)

        d = np.sqrt( dx_adj*dx_adj + dy_adj*dy_adj)
        R = 6378140 # m

        lat1 = cen_lat * (np.pi/180) # convert to radians
        lon1 = (cen_lon+180) * (np.pi/180) # convert to radians

        # Using haversine formula    
        bearing = (90 - math.degrees(math.atan2(dy_adj,dx_adj))) % 360
        bearing_rad = bearing * (np.pi/180)
        #print(bearing)

        lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(bearing_rad) )
        lon2 = lon1 + math.atan2( math.sin(bearing_rad)*math.sin(d/R)*math.cos(lat1), math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

        Nlat[i] = lat2 * (180/np.pi)
        Nlon[i] = lon2 * (180/np.pi) - 180

    return Nlat,Nlon

def plot_domain_boxes(ax,lat_p,lon_p,name,color):
    for i in np.arange(0,3):
        ax.plot([lon_p[name][i], lon_p[name][i+1]], [lat_p[name][i], lat_p[name][i+1]], color=color)
    ax.plot([lon_p[name][3], lon_p[name][0]], [lat_p[name][3], lat_p[name][0]], color=color)
    return()

# Calculate new domain lat and lon locations
lat_p[name],lon_p[name] = domain_box_from_center(cen_lat,cen_lon,ny*dy/2,nx*dx/2)

# Terrain levels
tlvls = np.arange(0,3001,500)

if plt_mult_dom == 0:

    # use glob to get list of h5 and header files
    h5filepath = path+'a-'+filetype+'*g'+gridnum+'.h5'
    h5files1 = sorted(glob.glob(h5filepath))
    hefilepath = path+'a-'+filetype+'*head.txt'
    hefiles1 = sorted(glob.glob(hefilepath))
    print(h5files1)
    
    # Loop through files
    for i in np.arange(0,len(h5files1)):   
        print(i)
    
        # Loop through variables
        for v in np.arange(0,len(vars)):
            var = vars[v]        
    
            # Loop through vertical levels to plot
            for k in np.arange(0,len(lvls)):
       
                znow = lvls[k] # get current altitude
                
                # Added plot function here for parallelization                
                plot_single_dom(var, znow, hefiles1[i], h5files1[i], x0,x1,y0,y1, scale, clvls, cmaps, tlvls, savepath, sname, gridnum)

elif plt_mult_dom == 1:

    # Find files for different grids and header files
    h5filepath = path+'a-A*g3.h5'
    h5files3 = sorted(glob.glob(h5filepath))
    
    h5filepath = path+'a-A*g2.h5'
    h5files2 = sorted(glob.glob(h5filepath))
    
    h5filepath = path+'a-A*g1.h5'
    h5files1 = sorted(glob.glob(h5filepath))
    
    hefilepath = path+'a-A*head.txt'
    hefiles1 = sorted(glob.glob(hefilepath))
    print(h5files1)
    
    for i in np.arange(0,len(h5files1)):
        print(i)
    
        # Loop through variables
        for v in np.arange(0,len(vars)):
            var = vars[v]
    
            # Loop through vertical levels to plot
            for k in np.arange(0,len(lvls)):
    
                znow = lvls[k] # get current altitude
    
                # Get RAMS vertical coodinate levels
                zcoords = np.array(rams_tools.calc_zcoords(hefiles1[i]))
                dzcoords = np.diff(zcoords)
                # Get the vertical coordinate index associated with the specified altitudes
                zid = np.where(np.abs(zcoords-znow) == np.min(np.abs(zcoords-znow)))[0][0]
                print(zid)
    
                # convert altitude to string for filenames
                lvl = int(np.array(zcoords[zid]))
                if lvl < 10:
                   lvlstr = '0000'+str(lvl)
                elif lvl < 100:
                   lvlstr = '000'+str(lvl)
                elif lvl < 1000:
                   lvlstr = '00'+str(lvl)
                elif lvl < 10000:
                   lvlstr = '0'+str(lvl)
                else:
                   lvlstr = str(lvl)
    
                # Grab current time from rams filename
                cur_time = os.path.split(h5files1[i])[1][9:21]
                print(cur_time)
    
                # Just open topography to get 2D dimensions for determining figure size
                rams_file = h5py.File(h5files1[0], 'r')
                topt = np.array(rams_file['TOPT'][:,:])
                fsx = 8;
                fsy = fsx*np.shape(topt)[0]/np.shape(topt)[1]/1.2
                # Make Figure
                fig,ax = plt.subplots(1,1,figsize=(fsx,fsy))
                rams_file.close()
    
                # loop through each potential grid
                for ii in np.arange(0,3):
    
                    # see if grid is present for current time, and if so, open up file
                    if ii == 0:
    
                       index = [idx for idx, s in enumerate(h5files1) if cur_time in s]
                       print(index)
                       rams_file = h5py.File(h5files1[index[0]], 'r')
                       gname = '1';
    
                    elif ii == 1:
                       if not len(h5files2):
                          continue
    
                       index = [idx for idx, s in enumerate(h5files2) if cur_time in s]
                       print(index)
                       if not len(index):
                          continue
    
                       rams_file = h5py.File(h5files2[index[0]], 'r')
                       gname = '2'
                    elif ii == 2:
                       if not len(h5files3):
                          continue
    
                       index = [idx for idx, s in enumerate(h5files3) if cur_time in s]
                       print(index)
                       if not len(index):
                          continue
    
                       rams_file = h5py.File(h5files3[index[0]], 'r')
                       gname = '3'
    
                     # Use zoomed in region; else use entire domain
                    if x0 >= 0:
                       lats = np.array(rams_file['GLAT'][x0:x1,y0:y1])
                       lons = np.array(rams_file['GLON'][x0:x1,y0:y1])
                       topt = np.array(rams_file['TOPT'][x0:x1,y0:y1])
                       if len(np.shape(rams_file[var])) == 2:
                           pltvar = np.array(rams_file[var][x0:x1,y0:y1])
                       else:
                           pltvar = np.array(rams_file[var][zid,x0:x1,y0:y1])
                    else:
                       lats = np.array(rams_file['GLAT'][:,:])
                       lons = np.array(rams_file['GLON'][:,:])
                       topt = np.array(rams_file['TOPT'][:,:])
                       if len(np.shape(rams_file[var])) == 2:
                           pltvar = np.array(rams_file[var][:,:])
                       else:
                           pltvar = np.array(rams_file[var][zid,:,:])
                   
                    # Multiple grids must use latlon dimensions; otherwise cannot overlap grids 
                    sname2 = 'LatLon'
    
                    a = ax.contourf(lons,lats,pltvar*scale[var],levels=clvls[var],cmap=cmaps[var],extend='both')
                    b = ax.contour(lons,lats,topt,levels=tlvls,colors='k')
    
                    # Plot outline of domains
                    ax.plot([np.nanmax(lons),np.nanmax(lons)],[np.nanmin(lats),np.nanmax(lats)],'-k')
                    ax.plot([np.nanmin(lons),np.nanmin(lons)],[np.nanmin(lats),np.nanmax(lats)],'-k')
                    ax.plot([np.nanmin(lons),np.nanmax(lons)],[np.nanmax(lats),np.nanmax(lats)],'-k')
                    ax.plot([np.nanmin(lons),np.nanmax(lons)],[np.nanmin(lats),np.nanmin(lats)],'-k')
    
                plt.colorbar(a,ax=ax)
    
                # Plot potential new domain
                if plt_new_dom == 1:
                    plot_domain_boxes(ax,lat_p,lon_p,name,'green')
                    sname3 = 'AddNew_'+str(cen_lat)+'_'+str(cen_lon)+'_'+str(ny)+'_'+str(nx)+'_'

                ax.set_title(var+' @ '+str(np.round(zcoords[zid]/1000,1))+' km:'+cur_time)
                ax.grid()
                plt.tight_layout()
    
                plt.savefig(savepath+sname2+sname+'_MultiGrid_'+sname3+var+'_'+lvlstr+'_'+cur_time+'.png')
                plt.close(fig)
