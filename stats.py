#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 09:22:56 2018

@author: therry

Module contenant les outils utilisés pour la création des fichiers de stats
"""

import gsw
import numpy as np
import netCDF_form as ncform
from scipy import interpolate as ip

import param as param
import general_tools as tools
import tile as ti


var_choice = {'zmean': ['NBbar', 'SAbar', 'CTbar', 'Ribar', 'BVF2bar'],
              'zstd': ['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std'], 
              'zdz': ['DZmean', 'DZstd', 'DZskew', 'EAPE']}


def create_stat_file(itile, atlas_infos):
    """Create statistics netcdf file
    
    :rtype: None
    """
    
    filename = '%s/stats_%003i.nc' % (param.get_path('stats'), itile)
    grid_lat, grid_lon = grid_coordinate(itile, atlas_infos['RESO'])
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)
    nlat, nlon = np.shape(lon_deg)

    print('create stat file %s' % filename)
    ncform.create_dim(filename, param.zref, nlat, nlon, atlas_infos)
    ncform.create_var(filename, var_choice[atlas_infos['TYPESTAT']])


def write_stat_file(itile, atlas_infos):
    """Write statistics into a netcdf file
    
    :rtype: None
    """

    filename = '%s/stats_%003i.nc' % (param.get_path('stats'), itile)
    res = {}

    res = compute_at_zref(itile, atlas_infos['RESO'], atlas_infos['MODE'], atlas_infos['DATE'], atlas_infos['TYPESTAT'])

    res['zref'] = param.zref
    ncform.write_var(filename, var_choice[atlas_infos['TYPESTAT']], res)

def read_stat_file(itile, typestat, reso, timeflag, date, mode, var_choice):
    """Read statistics into a netcdf file
    
    :rtype: dict"""

    filename = '%s/stats_%003i.nc' % (param.get_path('stats'), itile)
    print('read stat file : %s' % filename)

    res = ncform.read_var(filename, var_choice)
    return res


def date_mode_filter(mode, date, itile):
    """
    Make the tile filter to choose keep only the chosen mode ('R', 'A', 'D', 'AD' or 'RAD')
    and the profiles under the chosen date (year, month, day)
    Return the tile_extract according to the filters used.
    
    :rtype: dic
    """
    key_extraction = ['RHO', 'SA', 'CT', 'BVF2']
    zref_profiles = ti.read_zref_profiles(itile)
    argo_tile = ti.read_argo_tile(itile)
    julday = tools.conversion_gregd_juld(date)
    mode_list = list(mode)
    tile_extract = {}
    
    idx = []
    
    for m in mode_list:
        sub_argo_tile = argo_tile[(argo_tile['DATA_MODE'] == m) & (argo_tile['JULD'] < julday) & (argo_tile['FLAG'] == 0)]
        idx += list(sub_argo_tile.index)
    
    
    for k in key_extraction:
        tile_extract[k] = zref_profiles[k].loc[idx]

    return tile_extract


def compute_at_zref(itile, reso_deg, mode, date, typestat, tile_dict=None):
    """Compute the variables at depths zref
    
    :rtype: dict"""

    argo_tile = ti.read_argo_tile(itile)

    if tile_dict != None:
        tile = tile_dict
    else:
        tile = date_mode_filter(mode, date, itile)
    
    
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    # Trouver parade recherche d'indice NaN variables
    #  nanidx = np.where(np.isnan(CT) | np.isnan(SA))
    lat = argo_tile.loc[tile['CT'].index, 'LATITUDE']
    lon = argo_tile.loc[tile['CT'].index, 'LONGITUDE']
    
    
    nanidx = []
    for z in range(len(param.zref)):
        # On recherche les indices des lignes contenant des NaN classés par colonne
        # Chaque sous liste correspond à une colonne de la DataFrame
        nanidx.append(CT.index[CT.iloc[:, z].isnull()])

    grid_lat, grid_lon = grid_coordinate(itile, reso_deg)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)

    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(param.zref)
    nbprof = len(CT)    

    variables = {}

    # gridded arrays of CT, SA et RI means
    for i, v in enumerate(var_choice['zmean']):
        variables[v] = np.zeros((nz, nlat, nlon))

    for k in range(nbprof):
        #  print('%4i/%i' % (k, nbprof))
        # todo: weigh in time using juld,
        # e.g. only winter statistics
        time_weight = 1.
        xlon_rad = np.deg2rad(lon.iloc[k])
        xlat_rad = np.deg2rad(lat.iloc[k])
        weight = tools.compute_weight(lon_rad, lat_rad,
                                      xlon_rad, xlat_rad,
                                      reso_rad)
        weight *= time_weight
        for l in range(nz):
            if np.isnan(CT.iloc[k, l]) or np.isnan(SA.iloc[k, l]):
                pass
            else:
                variables['NBbar'][l, :, :] += weight
                variables['CTbar'][l, :, :] += weight*CT.iloc[k, l]
                variables['SAbar'][l, :, :] += weight*SA.iloc[k, l]
                variables['Ribar'][l, :, :] += weight*RI.iloc[k, l]
                variables['BVF2bar'][l, :, :] += weight*BVF2.iloc[k, l]

    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    coef = 1./variables['NBbar']
    coef[variables['NBbar'] < 1] = np.NaN

    variables['CTbar'] *= coef
    variables['SAbar'] *= coef
    variables['Ribar'] *= coef
    variables['BVF2bar'] *= coef


    if typestat =='zstd' or typestat == 'zdz':
        xlon_rad = np.deg2rad(lon)
        xlat_rad = np.deg2rad(lat)
        for i, v in enumerate(var_choice[typestat]):
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
                        weight = tools.compute_weight(lon_rad[j, i], lat_rad[j, i],
                                                      xlon_rad, xlat_rad,
                                                      reso_rad)
                        weight *= time_weight
                        drho = RI - variables['Ribar'][:, j, i]
                        dbvf2 = BVF2 - variables['BVF2bar'][:, j, i]
                        dCT = CT - variables['CTbar'][:, j, i]
                        interpolator = ip.interp1d(variables['Ribar'][:,j,i], param.zref, bounds_error=False)
                        p = gsw.p_from_z(-param.zref, lat[j])
                        g = gsw.grav(lat[j], p)
                        cs = gsw.sound_speed(variables['SAbar'][:, j, i], variables['CTbar'][:, j, i], p)
                        rho0 = variables['Ribar'][:, j, i].copy()
                        zrho = interpolator(RI)
                        dzstar = zrho-param.zref
                        dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
                        dSA = SA - variables['SAbar'][:, j, i]

                        weight = weight[:, np.newaxis] + np.zeros_like(param.zref)
                        weight[np.where(np.isnan(dz) | np.isnan(drho) | np.isnan(dCT) | np.isnan(dSA))] = 0.
                        #  weight[nanidx] = 0.
                        def average(field):
                            return np.nansum(weight*field, axis=0)
                        if typestat == 'zstd':
                            variables['CTstd'][:, j, i] = average(dCT**2)
                            variables['SAstd'][:, j, i] = average(dSA**2)
                            variables['BVF2std'][:, j, i] = average(dbvf2**2)
                            variables['Ristd'][:, j, i] = average(drho**2)
                               
                            coef = 1./(variables['NBstd']-1)
                            coef[variables['NBstd'] < 2] = np.nan

                                

                        if typestat == 'zdz':
            
                            variables['DZmean'][:, j, i] = average(dz)
                            variables['DZstd'][:, j, i] = average(dz**2)
                            variables['DZskew'][:, j, i] = average(dz**3)
                            variables['EAPE'][:, j, i] = average(dz*drho)
                                
                                
    if typestat == 'zstd':
        variables['CTstd'] = np.sqrt(coef*variables['CTstd'])
        variables['SAstd'] = np.sqrt(coef*variables['SAstd'])
        variables['Ristd'] = np.sqrt(coef*variables['Ristd'])
        variables['BVF2std'] = np.sqrt(coef*variables['BVF2std'])

    elif typestat == 'zdz':
        variables['DZmean'] *= coef
        variables['DZstd'] = np.sqrt(coef*variables['DZstd'])
        variables['DZskew'] *= coef/variables['DZstd']**3
        variables['EAPE'] *= 0.5*coef

    variables['lat'] = lat_deg
    variables['lon'] = lon_deg

    return variables


def grid_coordinate(itile, reso):
    """ 
    Returns the coordinates of each point of the grid for a given tile

    coordinates are round multiples of reso_deg
    reso sets the grid resolution, typically 0.5deg
    
    :rtype: numpy.ndarray, numpy.ndarray""" 

    tiles = ti.tile_definition()

    i = itile // tiles['NLON']
    j = itile % tiles['NLON']

    latmin = np.ceil(tiles['LATITUDE'][i]/reso)*reso
    latmax = np.floor(tiles['LATITUDE'][i+1]/reso)*reso
    lonmin = np.ceil(tiles['LONGITUDE'][j]/reso)*reso
    lonmax = np.floor(tiles['LONGITUDE'][j+1]/reso)*reso
    if (lonmax % 1. == 0.):
        # remove rightmost point (it's the leftmost point of the next
        # tile). That's also true for +180 (that is == -180).
        lonmax -= reso

    grid_lat = np.arange(latmin, (latmax+reso), reso)
    grid_lon = np.arange(lonmin, (lonmax+reso), reso)

    return(grid_lat, grid_lon)

def main(itile):
    """
    Main function of stats.py
    """
    create_stat_file(itile, param.atlas_infos)
    write_stat_file(itile, param.atlas_infos)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    atlas_infos = param.get_atlas_infos()
    create_stat_file(125, atlas_infos)
    write_stat_file(125, atlas_infos)