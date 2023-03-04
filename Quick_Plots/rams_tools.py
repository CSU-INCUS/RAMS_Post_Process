# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 08:34:21 2015

@author: sfreeman
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import h5py
import scipy.special

def make_ticks(minv, maxv, grid=False):
    div = (maxv-minv)//6
    l1 = [10, 7.5, 5, 2.5, 1]
    diff = math.pow(10,int(math.log(div, 10)))
    l2 = [diff * x for x in l1]
    td = min(l2, key=lambda x: abs(x-div))

    diff = minv//10
    l3 = [diff * x for x in [10,7.5, 5, 2.5, 0]]
    minv_new = min(l3, key=lambda x: abs(x-minv))
    diff = math.pow(10,int(math.log(maxv, 10)))       
    
    l3 = [diff * x for x in [10,7.5, 5, 2.5, 0]]
    #print(l3)

    maxv_new = min((i for i in l3 if (i-maxv)>0), key = lambda x: abs(x-maxv))
    
    # "safety valve" if large number, but low amount
    if 10*td < minv:
        minv_new = minv
        maxv_new = maxv
    
    if grid is True:
         if round(td)!=td:
                td = max(10,1, key = lambda x: x-td)
    
    
    #print(minv_new,maxv_new, td)
    return np.arange(minv_new, maxv_new+1, td)

def make_rams_plot(cfp = None, csp = None, quiverx = None, quivery = None, skip=1, 
                           minx=0, maxx=420, miny=0, maxy=510, 
                           zlev=0, tnum=0, title='', clevs = None, cslevs = None, cmap = 'jet',
                           axis = None, display_colorbar = True, showgrid = False,
                           extend_white=False, ticks = None, streamlines = False, 
                           xticks = None, yticks = None, cbar_title='', f=None, cscols='k', minutes_only=False, cslw=2):
    '''
    Default plotter for RAMS input data.
    
    Args:
        cfp: a 2-D array of values to plot a filled contour
        csp: a 2-D array of values to plot a non-filled contour
        quiverx: the x values of a quiver plot
        quivery: the y values of a quiver plot
        skip: the amount to skip by in each direction in a quiver plot
        minx: the minimum x-coordinate you are plotting. Note that 
              this *must* match up with the cfp and csp given. 
              this function will not subset your data automatically!
        maxx: the maximum x-coordinate
        miny: the minimum y-coordinate
        maxy: the maximum y-coordinate
    '''


    x = np.arange(minx,maxx)
    y = np.arange(miny,maxy)
    X,Y = np.meshgrid(x,y)
    if axis is not None:
        plotting = axis
    else:
        plotting = plt.axes()

    plotting.set_aspect('equal')
    plotting.set_xticks(make_ticks(minx, maxx, grid=True))
    plotting.set_yticks(make_ticks(miny, maxy, grid=True))
    if xticks is not None and yticks is not None:
        plotting.set_xticks(xticks)
        plotting.set_yticks(yticks)
    if xticks is not None:
        plotting.set_xticks(xticks)


            
    # filled contour plots
    if cfp is not None:
        cs = plotting.contourf(X, Y, cfp, extend='both', cmap = plt.get_cmap(cmap), levels = clevs)
        if display_colorbar is not False:
            div = make_axes_locatable(plotting)
            cax = div.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(cs, cax=cax, format="%.4g")
            if ticks is not None:
                cbar.set_ticks(ticks)
            if cbar_title !='':
                cbar.set_label(cbar_title)
        if extend_white is True:
            cs.cmap.set_under('white')

    if csp is not None:
        try:
            cnf = plotting.contour(X, Y, csp, 40, colors=cscols, levels=cslevs, lw = cslw)
            plotting.clabel(cnf, inline=1, fontsize=10) 
        except ValueError:
            print("No contours on this plot due to no variability")
    
    if quiverx is not None and quivery is not None:
        if streamlines is True:
            plotting.streamplot(X[::skip,::skip], Y[::skip,::skip], 
               quiverx,quivery, splotdensity=2, color='k')
        else:
            plotting.quiver(X[::skip,::skip], Y[::skip,::skip],
                        quiverx[::skip,::skip],quivery[::skip,::skip])

    if showgrid == True:
        plotting.grid()
    
    if minutes_only==False:
        plotting.set_title(title+"\nat z = "+ str(zlev) + 
                  '  and t = '+str(tnum)+ ' minutes')
    else:
        plotting.set_title(str(tnum)+" minutes")
    return cs

def calc_vort(us, vs, dx, dy):
    '''
    Calculates relative vorticity in the z direction 
    us: u winds
    vs: v winds
    dx: grid spacing (in m) in x direction
    dy: grid spacing (in m) in y direction
    '''
    dudy, dudx = np.gradient(us)
    dvdy, dvdx = np.gradient(vs)
    dudy = dudy/dy
    dvdx = dvdx/dx

    return dvdx - dudy

def calc_div(us, vs, dx, dy):
    '''
    Calculates divergence  
    us: u winds
    vs: v winds
    dx: grid spacing (in m) in x direction
    dy: grid spacing (in m) in y direction
    '''
    dudy, dudx = np.gradient(us)
    dvdy, dvdx = np.gradient(vs)
    dudx = dudx/dx
    dvdy = dvdy/dy

    return dudx+dvdy

def calc_wsp(us, vs):
    '''
    Calculates wind speed.
    us: u winds
    vs: v winds
    '''
    return np.power(np.power(us,2)+np.power(vs,2), 0.5)

def calc_zcoords(fname):
    '''
    Calculates/reads in z coordinates from RAMS file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__ztn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords

def calc_zmn(fname):
    '''
    Calculates/reads in z coordinates from RAMS file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__zmn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords

def calc_ztn(fname):
    '''
    Calculates/reads in ztn coordinates from RAMS file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__ztn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords


def calc_pi01dn01(fname):
    '''
    Reads base state pi from the RAMS header file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__pi01dn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords

def calc_th01dn01(fname):
    '''
    Reads base state theta from the RAMS header file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__th01dn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords

def calc_dzmn01(fname):
    '''
    Reads base state theta from the RAMS header file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__dzmn01')
    numlines = int(mylist[ix+1])
    zcoords = mylist[ix+2:ix+2+numlines]
    zcoords = [float(x) for x in zcoords]
    return zcoords


def calc_xcoords(fname):
    '''
    Calculates/reads in x coordinates from RAMS file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__xtn01')
    numlines = int(mylist[ix+1])
    xcoords = mylist[ix+2:ix+2+numlines]
    xcoords = [float(x) for x in xcoords]
    return xcoords

def calc_ycoords(fname):
    '''
    Calculates/reads in y coordinates from RAMS file.
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__ytn01')
    numlines = int(mylist[ix+1])
    ycoords = mylist[ix+2:ix+2+numlines]
    ycoords = [float(x) for x in ycoords]
    return ycoords


def calc_bstden(fname):
    '''
    Calculates/reads in the base state density
    '''
    with open(fname) as f:
        mylist = f.read().splitlines() 
    ix = mylist.index('__dn01dn01')
    numlines = int(mylist[ix+1])
    bstden = mylist[ix+2:ix+2+numlines]
    bstden = [float(x) for x in bstden]
    return bstden

def get_file_list(fnames):
    '''
    Gets a list of RAMS files as a list of h5py files.
    fnames: a list of file names (likely from glob.glob)
    '''
    f_ret = list()
    for fname in sorted(fnames):
        f_ret.append(h5py.File(fname, 'r'))
    return f_ret

def calc_watersat(xpress, xtemp):
    '''
    Returns the water saturation vapor pressure for a temp and pressure.
    '''
    raise NotImplementedError


def calc_relh(rv, pi, theta):
    '''
    Calculates the relative humidity given RAMS model output
    a: RV
    b: PI 
    c: THETA
    '''
    raise NotImplementedError
    cp = 1004
    cpor = 00
    p00 = 00
    xtemp=theta*pi/cp
    xpress=(pi/cp)**cpor*p00
    return 100.*rv/rsatmix(xpress,xtemp)

def dropsizedist(diameter, meandiameter, nu):
    '''
    Calculates the drop size distribution function value at a given diameter
    given the shape function and mean diameter.
    '''

    dn = meandiameter / (np.power(scipy.special.gamma(nu+3)/scipy.special.gamma(nu),1./3))
    #print(dn, nu, meandiameter)
    top = (np.power(diameter,(nu-1))*np.exp(-1*diameter/dn))
    bottom = scipy.special.gamma(nu)*np.power(dn,nu)
    #print(top,bottom)
    return top/bottom



def plot_plumes(axes, files, rvars, runnames, vartype_names):
    for i, var in enumerate(rvars):
        axes[i].plot(files)

def calc_thetarho(RV, RTP, theta):
    return theta * ((1+RV/0.622)/(1+RTP))

def calc_coldpool_strength(RV, RVO, RTP, theta, theta_mean, thetarho_mean, dzarray, zt_cond_B):
    '''
    Calculates c^2
    RV: RV from model (z, y, x array)
    RVO: RV Naught, initial RV
    RTP: RTP from model (z, y, x array)
    theta: theta from model (z, y, x array)
    thetarho_mean: Mean theta rho array
    dzarray: dz array
    zt_cond_B: the top of the buoyancy condition
    '''
    thetarho_prime = calc_thetarho(RV, RTP, theta) - thetarho_mean
    tho = theta_mean * (1 + 0.61 * rvo)
    buoyancy = 9.81 * (())

    return -2* 9.81 * np.ma.sum((thetarho_prime/thetarho_mean) * dzarray, axis=0)

def calc_press_from_pi(pi):
    '''
    Calculates pressure from the Pi parameter in RAMS
    
    Inputs:
    pi: pi parameter from RAMS
    
    Outputs:
    pressure in hPa
    '''
    #divide by cp
    pi = pi/1004
    
    pi = np.power(pi, 1004/287)
    return pi*1000

def calc_press_from_pi_dask(pi):
    import dask as da
    '''
    Calculates pressure from the Pi parameter in RAMS
    
    Inputs:
    pi: pi parameter from RAMS
    
    Outputs:
    pressure in hPa
    '''
    #divide by cp
    pi = pi/1004
    
    pi = pi**(1004/287)
    return pi*1000


def calc_T_from_theta(theta, pi):
    '''
    Calculates air temperature from theta and pi in RAMS
    
    Inputs:
    theta: potential temperature (K)
    pi: exner function (RAMS; J/kg*K)
    
    Outputs:
    air temperature in K
    '''
    
    return (pi/1004)*theta


def calc_vapor_press(pi, rv):
    '''
    Calculates the vapor pressure using the Hobbs 1977 (pg 71) formula
    Inputs:
    pi: exner function (RAMS; J/kg*K)
    rv: water vapor mixing ratio (kg/kg)
    
    Outputs:
    the vapor pressure (e) in hPa
    '''
    epsilon = 18.01528/28.9644
    press = calc_press_from_pi(pi)
    
    return press * rv/(epsilon+rv)

def saturation_vapor_pressure(theta, pi):
    '''
    Calculates the saturation vapor pressure using the Bolton1980 formula converted to K by metpy
    Inputs:
    theta: potential temperature (K)
    pi: exner function (RAMS; J/kg*K)

    Outputs: 
    the saturation vapor pressure (es) in hPa
    '''
    temp = calc_T_from_theta(theta, pi)
    return 6.112 * np.exp(17.67 * (temp - 273.15)
                                    / (temp - 29.65))

    

def calc_dewpt_from_theta_rv_dask(theta, rv, pi):
    import dask as da
    '''
    Calculates dewpoint from theta, rv, and pi in RAMS
    Inputs:
    theta: potential temperature (K)
    pi: exner function (RAMS; J/kg*K)
    rv: water vapor mixing ratio (kg/kg)
    '''
    
    e = calc_vapor_press(pi, rv)
    
    val =da.array.log(e /  6.112)
    return 243.5 * val / (17.67 - val)


def prep_rams_iris_cube_tobac_3D(incube, xcoords, ycoords):
    '''
    Preps the input rams cube converted from xarray (incube) and get it ready for tobac.
    
    Inputs:
    incube: the input cube (usually of one variable)
    xcoords: x-coordinates from xtn
    ycoords: y-coordinates from ytn
    
    Outputs:
    none, works inplace. 
    '''
    from iris import coords
    incube.coords()[0].rename('time')
    
    coord_system=None
    x_coord=coords.DimCoord(np.arange(len(xcoords)), long_name='x', units='1', bounds=None, attributes=None, coord_system=coord_system)
    incube.add_dim_coord(x_coord,3)
    y_coord=coords.DimCoord(np.arange(len(ycoords)), long_name='y', units='1', bounds=None, attributes=None, coord_system=coord_system)
    incube.add_dim_coord(y_coord,2)
    
    projection_x_coord=coords.DimCoord(xcoords, standard_name='projection_x_coordinate', long_name='x', var_name='x', units='m', bounds=None, attributes=None, coord_system=coord_system)
    incube.add_aux_coord(projection_x_coord,(3))
    projection_y_coord=coords.DimCoord(ycoords, standard_name='projection_y_coordinate', long_name='y', var_name='y', units='m', bounds=None, attributes=None, coord_system=coord_system)
    incube.add_aux_coord(projection_y_coord,(2))

