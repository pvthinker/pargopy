# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 09:34:40 2018

@author: herry
"""
from netCDF4 import Dataset
from scipy import interpolate as ip
import os
import numpy as np
import gsw as gsw
import time
import pickle as pickle
import general_tools as tools
import research_tools as research
import argodb as argo
import argotools as argotools
import interpolation_tools as interpolation
import tile as tiler


zref = argotools.zref
path_localdata = '/local/tmp/1/herry/pargopy/stats/'


#  ----------------------------------------------------------------------------
def creating_stats(itile, typestat, reso_deg, timeflag):
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    if typestat == 'zmean':
            res = compute_mean_at_zref(itile, reso_deg)
            write_stats(res, itile, typestat, reso_deg, timeflag)
    if typestat == 'zstd':
            res = compute_std_at_zref(itile, reso_deg)
            write_stats(res, itile, typestat, reso_deg, timeflag)


#  ----------------------------------------------------------------------------
def write_stats(stats, itile, typestat, reso_deg, timeflag):
    """Write statistics into a netcdf file"""
    # idem, depend du type de stat
    filename = '/local/tmp/1/herry/pargopy/stats/%s_%s_%s_%s.pkl' % (typestat, reso_deg, timeflag, itile)
    with open(filename, 'w') as f:
        pickle.dump(stats, f)


#  ----------------------------------------------------------------------------
def read_stats(itile, typestat, reso_deg, timeflag):
    filename = '/local/tmp/1/herry/pargopy/stats/%s_%s_%s_%s.pkl' % (typestat, reso_deg, timeflag, itile)
    print(filename)
    with open(filename, 'r') as f:
        stats = pickle.load(f)
    return stats


#  ----------------------------------------------------------------------------
def compute_mean_at_zref(itile, reso_deg):
    """Compute the mean at depths zref"""
    tile = tiler.read_tile(itile)
    argodic = research.read_argo_filter(itile)
    argodb = argo.read_argodb()
    output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    CT, SA, RI, lat, lon, juld, zref = tile['CT'], tile['SA'], tile['RHO'], tile['LATITUDE'], tile['LONGITUDE'], tile['JULD'], tile['ZREF']
    nbprof = output['IPROF']
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso_deg)

    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(zref)

    # gridded arrays of CT, SA et RI means
    NBbar = np.zeros((nz, nlat, nlon))
    CTbar = np.zeros((nz, nlat, nlon))
    SAbar = np.zeros((nz, nlat, nlon))
    RIbar = np.zeros((nz, nlat, nlon))
    for k in range(len(nbprof)):
        # todo: weigh in time using juld,
        # e.g. only winter statistics
        time_weight = 1.

        xlon_rad = np.deg2rad(lon[k])
        xlat_rad = np.deg2rad(lat[k])
        weight = tools.compute_weight(lon_rad, lat_rad,
                                      xlon_rad, xlat_rad,
                                      reso_rad)
        weight *= time_weight

        for l in range(nz):
            if np.isnan(CT[k, l]) or np.isnan(SA[k, l]):
                pass
            else:
                NBbar[l, :, :] += weight
                CTbar[l, :, :] += weight*CT[k, l]
                SAbar[l, :, :] += weight*SA[k, l]
                RIbar[l, :, :] += weight*RI[k, l]

    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    coef = 1./NBbar
    coef[NBbar < 1] = np.NaN

    CTbar *= coef
    SAbar *= coef
    RIbar *= coef

    res = {'lon_deg': lon_deg,
           'lat_deg': lat_deg,
           'NBbar': NBbar,
           'CTbar': CTbar,
           'SAbar': SAbar,
           'RIbar': RIbar}

    return res


def compute_std_at_zref(itile, reso_deg, timeflag):
    """Compute the standard deviations at depths zref"""

    # gridded arrays of CT, SA variances
    lon_deg, lat_deg, NBbar, CAbar, SAbar, RHObar = read_stats('zmean', itile, reso_deg, timeflag) # read it from the file
    tile = tiler.read_tile(itile)
    argodic = research.read_argo_filter(itile)
    argodb = argo.read_argodb()
    output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    CT, SA, RI, lat, lon, juld, zref = tile['CT'], tile['SA'], tile['RHO'], tile['LATITUDE'], tile['LONGITUDE'], tile['JULD'], tile['ZREF']
    nbprof = output['IPROF']
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(zref)
    
    xlon_rad = np.deg2rad(lon)
    xlat_rad = np.deg2rad(lat)

    CTstd = np.zeros((nz, nlat, nlon))
    SAstd = np.zeros((nz, nlat, nlon))
    DZstd = np.zeros((nz, nlat, nlon))
    DRHOstd = np.zeros((nz, nlat, nlon))
    EAPE = np.zeros((nz, nlat, nlon))

    wmin = 5e-3 # minimum weight below which a profile is drop
    
    # double loop on each grid point (instead of a loop on each profile)
    for j in range(nlat):
        for i in range(nlon):
            
            time_weight = 1.

            # do all profiles
            weight = tools.compute_weight(lon_rad[j, i], lat_rad[j, i],
                                            xlon_rad, xlat_rad,
                                            reso_rad)
            weight *= time_weight
            print(np.shape(RHObar))
            interpolator = ip.interp1d(RHObar, zref)

            p = gsw.p_from_z(-zref, lat[j, i])
            g = gsw.grav(lat[j, i], p)
            cs = gsw.sound_speed(SAbar[:, j, i], CAbar[:, j, i], p)
            # the copy is to make data contiguous in memory
            rho0 = RHObar[:, j, i].copy()
        
        for l in range(nbprof):
            if weight[l] < wmin:
                pass
            else:
                # the vertical isopycnal displacement dz
                # is computed by first
                # interpolating RI(z) on the reference profile
                # z0(rho) [inverse mapping of rho0(z)]
                # we get dzstar
                # then by adding a correction due to
                # compressibility
                dzstar = interpolator(RI[l, :]) - zref
                drho = RI[l, :] - rho0
                # add correction
                dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
                eape = 0.5*dz*drho
                
                for k in range(nz):
                    if np.isnan(dz[k]) or np.isnan(drho[k]):
                        pass
                    else:
                        CTstd[k, j, i] += weight*CT[l, k]**2
                        SAstd[k, j, i] += weight*SA[l, k]**2
                        DZstd[k, j, i] += weight*dz[k]**2
                        DRHOstd[k, j, i] += weight*drho[k]**2
                        EAPE[k, j, i] += weight*dz[k]*drho[k]
    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    # std = sqrt( 1/(n-1) sum_i (x_i -xbar)^2)
    #     = sqrt( 1/(n-1) sum_i x_i^2 - n/(n-1)*xbar^2)
    coef = 1./(NBbar-1)
    coef[NBbar < 2] = 0.

    CTstd = np.sqrt( coef*(CTstd-NBbar*CAbar)) 
    SAstd = np.sqrt( coef*(SAstd-NBbar*SAbar)) 
    DZstd = np.sqrt( coef*DZstd)
    DRHOstd = np.sqrt( coef*DRHOstd)
    EAPE = 0.5*coef*EAPE
    
    res = {'lon_deg': lon_deg,
           'lat_deg': lat_deg,
           'CTstd': CTstd,
           'SAstd': SAstd,
           'DZstd': DZstd,
           'DRHOstd': DRHOstd,
           'EAPE': EAPE}

    return res


def define_grid(minlon, maxlon, minlat, maxlat, reso_deg):
    """ setup the grid coordinates (in degrees)
    coordinates are round multiples of reso_deg
    reso_deg sets the grid resolution, typically 0.5deg""" 

    minlon = np.ceil(minlon/reso_deg)*reso_deg
    maxlon = np.floor(maxlon/reso_deg)*reso_deg

    minlat = np.ceil(minlat/reso_deg)*reso_deg
    maxlat = np.floor(maxlat/reso_deg)*reso_deg

    lon1D = np.arange(minlon, maxlon+reso_deg, reso_deg)
    lat1D = np.arange(minlat, maxlat+reso_deg, reso_deg)

    lon, lat = np.meshgrid(lon1D, lat1D)

    return lon, lat


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    creating_stats(0, 'zmean', 0.5, 'annual')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)