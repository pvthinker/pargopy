#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 07:11:16 2018

@author: therry

Module contenant la série d'outils pour les tâches maitre/esclaves telles que :
    - Génération d'atlas
        - Interpolation à l'aide des outils d'interpolation (interpolation_tools.py)
        - Génération des statistiques (zmean, zstd, zdz, ...)
        - Propagation inverse des flag après traitement des données Argo
        - Collage des "dalles" de l'atlas

Ce module contient également les outils permettant la génération d'atlas sur
un domaine restreint :
    - Génération sous-domaine d'atlas avec 4 points du globe (longitude, latitude)
    - Génération sous-domaine dans une zone géographique définie ( Mer Méditerranée,
      Pacifique-Sud, ...)



"""
import os
from netCDF4 import Dataset
import numpy as np

import param as param
import tile as ti
import stats as stats
import netCDF_form as ncform
import general_tools as tools

atlas_infos = param.get_atlas_infos()
tiles = ti.tile_definition()

def gridindex2lonlat(ix, iy):
    
    lonmin_glo = -180.
    latmin_glo = -80.

    lon = lonmin_glo + ix*atlas_infos['RESO']
    lat = latmin_glo + iy*atlas_infos['RESO']
    return lon, lat

def get_glo_grid():
    lonsize = []
    longlo = []
    latsize = []
    latglo = []
    for i in range(tiles['NLON']):
        itile = i
        lat, lon = stats.grid_coordinate(itile, atlas_infos['RESO'])
        lonsize.append(len(lon))
        longlo += list(lon)

    for j in range(tiles['NLAT']):
        itile = j*tiles['NLON']
        lat, lon = stats.grid_coordinate(itile, atlas_infos['RESO'])
        latsize.append(len(lat))
        latglo += list(lat)

    lonsize = np.asarray(lonsize)
    longlo = np.asarray(longlo)
    latsize = np.asarray(latsize)
    latglo = np.asarray(latglo)
    return lonsize, latsize, longlo, latglo


def glue_tiles():
    """ Glue the stats tiles together into a global 3D atlas"""


    lonsize, latsize, lon, lat = get_glo_grid()

    zref = param.zref
    nz = len(zref)

    # get the list of variables of this atlas
    listvar = stats.var_choice[atlas_infos['TYPESTAT']]

    print('create atlas for variables: '+', '.join(listvar))

    ncfile = param.get_atlas_filename()
    ncform.create_dim(ncfile, zref, len(lat), len(lon), atlas_infos)
    ncform.create_var(ncfile, listvar)
    ncform.create_var(ncfile, ['zref', 'lon', 'lat'])

    lon2d, lat2d = np.meshgrid(lon, lat)

    nc = Dataset(ncfile, 'r+')
    nc.variables['zref'][:] = zref
    nc.variables['lon'][:, :] = lon2d
    nc.variables['lat'][:, :] = lat2d
    
    j0 = 0
    for j in range(tiles['NLAT']):
        i0 = 0
        for i in range(tiles['NLON']):
            itile = tools.ij2tile(i, j)
            i1 = i0+lonsize[i]
            j1 = j0+latsize[j]
            if itile % 30 == 0:
                print('glue tile %3i: lon=%7.2f - lat=%7.2f' % (itile, lon[i0], lat[j0]))
            
            fname = stats.generate_filename(itile)
            if os.path.isfile(fname):
                with Dataset(fname) as ncf:
                    for v in listvar:
                        z3d = ncf.variables[v][:, :, :]
                        for k in range(nz):
                            nc.variables[v][k, j0:j1, i0:i1] = z3d[k, :, :]
            else:
                for v in listvar:
                    for k in range(nz):
                        nc.variables[v][k, j0:j1, i0:i1] = 0.

            i0 = i1
        j0 = j1
    nc.close()
    print('Atlas %s has been generated successfully' % ncfile)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':

    glue_tiles()