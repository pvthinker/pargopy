#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 12:57:38 2018

@author: therry
"""

# File used to know how to calculate each variable of stats

import numpy as np
import argotools as tools
import stats as stats
import general_tools as general


zref = tools.zref

def compute_at_zref(itile, reso_deg, mode, date, var_choice):
    """Compute the variables at depths zref
    
    :rtype: dict"""
    tile = stats.date_mode_filter(mode, date, itile)
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
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
    typestat = []

    for i, v in enumerate(var_choice):
        variables[v] = np.zeros((nz, nlat, nlon))
        if v.find('bar') != -1:
            typestat.append('zmean')
        else:
            typestat.append('zstd')

        if typestat[i] == 'zmean':
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
                        if v != 'NBbar':
                            # NBbar is calculated before and in another way
                            variables[v][l, :, :] += weight*CT[k, l]
                        else:
                            pass
                coef = 1./variables['NBbar']
                coef[variables['NBbar'] < 1] = np.NaN
                if v != 'NBbar':
                    # NBbar is calculated before and in another way
                    variables[v] *= coef
                else:
                    pass
    print(variables)
    print(typestat)
    return variables

# Still the same to do with the std variables

#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    compute_at_zref(0, 0.5, 'D', ['2017', '12', '31'], ['NBbar', 'DZmean', 'CTbar', 'CTstd'])
