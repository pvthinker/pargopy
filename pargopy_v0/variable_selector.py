#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 12:57:38 2018

@author: therry
"""

# File used to know how to calculate each variable of stats

import numpy as np
from math import sqrt
import gsw as gsw
from scipy import interpolate as ip
import argotools as tools
import stats as stats
import general_tools as general


zref = tools.zref
block_choice = ['zmean', 'zstd', 'zdz']  # , 'rmean', 'rstd', 'rdz'
var_choice = {'zmean':['NBbar', 'SAbar', 'CTbar', 'Ribar', 'BVF2bar'],
              'zstd':['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std'], 
              'zdz':['DZmean', 'DZstd', 'DZskew', 'EAPE']}

def compute_at_zref(itile, reso_deg, mode, date, block_choice, tile_dict=None):
    """Compute the variables at depths zref
    
    :rtype: dict"""
    if tile_dict != None:
        tile= tile_dict
    else:
        tile = stats.date_mode_filter(mode, date, itile)
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    nanidx = np.where(np.isnan(CT) | np.isnan(SA))
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']
    grid_lat, grid_lon = stats.grid_coordinate(itile, reso_deg)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)

    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(zref)
    nbprof = len(CT)    

    variables = {}

    for b in block_choice:
        # gridded arrays of CT, SA et RI means
        for i, v in enumerate(var_choice['zmean']):
            variables[v] = np.zeros((nz, nlat, nlon))

        for k in range(nbprof):
            #  print('%4i/%i' % (k, nbprof))
            # todo: weigh in time using juld,
            # e.g. only winter statistics
            time_weight = 1.
            xlon_rad = np.deg2rad(lon[k])
            xlat_rad = np.deg2rad(lat[k])
            weight = general.compute_weight(lon_rad, lat_rad,
                                            xlon_rad, xlat_rad,
                                            reso_rad)
            weight *= time_weight
            for l in range(nz):
                if np.isnan(CT[k, l]) or np.isnan(SA[k, l]):
                    pass
                else:
                    variables['NBbar'][l, :, :] += weight
                    variables['CTbar'][l, :, :] += weight*CT[k, l]
                    variables['SAbar'][l, :, :] += weight*SA[k, l]
                    variables['Ribar'][l, :, :] += weight*RI[k, l]
                    variables['BVF2bar'][l, :, :] += weight*BVF2[k, l]

        # normalize with the number of profiles (fractional
        # because NBbar is fractionnal)
        coef = 1./variables['NBbar']
        coef[variables['NBbar'] < 1] = np.NaN

        variables['CTbar'] *= coef
        variables['SAbar'] *= coef
        variables['Ribar'] *= coef
        variables['BVF2bar'] *= coef


        if b =='zstd' or b == 'zdz':
            xlon_rad = np.deg2rad(lon)
            xlat_rad = np.deg2rad(lat)
            for i, v in enumerate(var_choice[b]):
                variables[v] = np.zeros((nz, nlat, nlon))
            variables['NBstd'] = variables['NBbar']

            if len(lat) == 0:
                pass
            else:
                for j in range(nlat):
                    for i in range(nlon):
                        if len(lat) < j+1:
                            pass
                        else:
                            time_weight = 1.
                            weight = general.compute_weight(lon_rad[j, i], lat_rad[j, i],
                                                            xlon_rad, xlat_rad,
                                                            reso_rad)
                            weight *= time_weight
                            drho = RI - variables['Ribar'][:, j, i]
                            dbvf2 = BVF2 - variables['BVF2bar'][:, j, i]
                            dCT = CT - variables['CTbar'][:, j, i]
                            interpolator = ip.interp1d(variables['Ribar'][:,j,i], zref, bounds_error=False)
                            p = gsw.p_from_z(-zref, lat[j])
                            g = gsw.grav(lat[j], p)
                            cs = gsw.sound_speed(variables['SAbar'][:, j, i], variables['CTbar'][:, j, i], p)
                            rho0 = variables['Ribar'][:, j, i].copy()
                            zrho = interpolator(RI)
                            dzstar = zrho-zref
                            dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
                            dSA = SA - variables['SAbar'][:, j, i]

                            weight = weight[:, np.newaxis] + np.zeros_like(zref)
                            weight[np.where(np.isnan(dz) | np.isnan(drho) | np.isnan(dCT) | np.isnan(dSA))] = 0.
                            weight[nanidx] = 0.
                            def average(field):
                                return np.nansum(weight*field, axis=0)
                            if b == 'zstd':
                                variables['CTstd'][:, j, i] = average(dCT**2)
                                variables['SAstd'][:, j, i] = average(dSA**2)
                                variables['BVF2std'][:, j, i] = average(dbvf2**2)
                                variables['Ristd'][:, j, i] = average(drho**2)
                               
                                coef = 1./(variables['NBstd']-1)
                                coef[variables['NBstd'] < 2] = np.nan

                                

                            if b == 'zdz':
            
                                variables['DZmean'][:, j, i] = average(dz)
                                variables['DZstd'][:, j, i] = average(dz**2)
                                variables['DZskew'][:, j, i] = average(dz**3)
                                variables['EAPE'][:, j, i] = average(dz*drho)
                                
                                
        if b == 'zstd':
            variables['CTstd'] = np.sqrt(coef*variables['CTstd'])
            variables['SAstd'] = np.sqrt(coef*variables['SAstd'])
            variables['Ristd'] = np.sqrt(coef*variables['Ristd'])
            variables['BVF2std'] = np.sqrt(coef*variables['BVF2std'])
            
        elif b == 'zdz':
            variables['DZmean'] *= coef
            variables['DZstd'] = np.sqrt(coef*variables['DZstd'])
            variables['DZskew'] *= coef/variables['DZstd']**3
            variables['EAPE'] *= 0.5*coef


    variables['lat'] = lat_deg
    variables['lon'] = lon_deg
    print(variables['CTstd'].min())
    print(variables['CTstd'].max())
    print(variables['SAstd'].min())
    print(variables['SAstd'].max())

    return variables

#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    compute_at_zref(52, 0.5, 'D', ['2017', '12', '31'], ['zstd'])
