# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:24:20 2018

@author: herry

Tools used to generate the filters that are used after to generate the tiles

"""
import os
import time
import numpy as np
import time
import pickle
import matplotlib.pyplot as plt
#  from mpl_toolkits.basemap import Basemap
import param as param
import argotools as argotools
import melted_functions as melted

path_to_filter = param.path_to_filter
path_to_data = param.path_to_data


#  ----------------------------------------------------------------------------
def tile_definition():
    """Creating the variables
    
    :rtype: float, float, int, int, float, float"""

    minlon = -180.
    maxlon = 180.
    latmin = -80.
    maxlat = 80.

    nlon = 20
    nlat = 15

    minmargin = 1.

    deltalon = maxlon-minlon
    deltalat = maxlat-latmin

    lon = minlon + np.arange(nlon+1) * deltalon/nlon
    lat = latmin + np.arange(nlat+1) * deltalat/nlat

    # distribute the latitudes so that their differences
    # vary in cos(lat)
    # do it via an iterative method
    for k in range(5):
        latm = 0.5*(lat[1:]+lat[:-1])
        dlat = np.diff(lat) * np.cos(latm*np.pi/180)
        dlat = dlat*0 + np.mean(dlat)
        dlat = dlat / np.cos(latm*np.pi/180)
        dlat = dlat/sum(dlat)*deltalat
        lat[1:] = latmin + np.cumsum(dlat)

    marginlat = minmargin / np.cos(latm*np.pi/180)
    marginlon = 2

    return lat, lon, nlat, nlon, marginlat, marginlon


#  ----------------------------------------------------------------------------
def creating_tiles():
    """Giving values to the variables
    
    :rtype: None"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    
    #  argodb = argotools.read_argodb()
    argodb = melted.read_dic('argodb', path_to_data)
    lat, lon, nlat, nlon, marginlat, marginlon = tile_definition()
    k = 0
    for i in range(nlat):
        for j in range(nlon):
            latmin = lat[i] - marginlat[i]
            latmax = lat[i + 1] + marginlat[i]
            lonmin = lon[j] - marginlon
            lonmax = lon[j + 1] + marginlon
            if lonmin < -180:
                lonmin += 360
            elif lonmax > 180:
                lonmax -= 360
            else:
                pass
            #  WITH_M : with margin
            #  NO_M : without margin
            res = {'LATMIN_WITH_M': latmin,
                   'LATMAX_WITH_M': latmax,
                   'LONMIN_WITH_M': lonmin,
                   'LONMAX_WITH_M': lonmax,
                   'LATMIN_NO_M': lat[i],
                   'LATMAX_NO_M': lat[i+1],
                   'LONMIN_NO_M': lon[j],
                   'LONMAX_NO_M': lon[j+1],
                   'MARGINLAT': marginlat[i],
                   'MARGINLON': marginlon}

            argo_extract = get_idx_from_tiles_lim(res, argodb)
            test_tiles(argo_extract, k)
            melted.write_dic('argodic%003i' % k, argo_extract, path_to_filter)
            #  write_argo_filter(argo_extract, k)
            k += 1


#  ----------------------------------------------------------------------------
def test_tiles(argo_extract, i):
    """ Test to know if the tiles are correctly done with the lat and lon 
    limits
    
    :rtype: None"""
    idx1 = np.where(argo_extract['LATITUDE'] > argo_extract['LATMAX_NO_M'] + argo_extract['MARGINLAT'])
    idx2 = np.where(argo_extract['LATITUDE'] < argo_extract['LATMIN_NO_M'] - argo_extract['MARGINLAT'])
    idx3 = np.where(argo_extract['LONGITUDE'] > argo_extract['LONMAX_NO_M'] + argo_extract['MARGINLON'])
    idx4 = np.where(argo_extract['LONGITUDE'] < argo_extract['LONMIN_NO_M'] - argo_extract['MARGINLON'])
    if (idx1[0] != []) | (idx2[0] != []) | (idx3[0] != []) | (idx4[0] != []):
        print('There is an error with the dimensions of the tile number %i' % i)
        exit(0)


#  ----------------------------------------------------------------------------
def write_argo_filter(argo_extract, i):
    """Writing each filter corresponding to the future tiles
    
    :rtype: None"""
    with open('%s/argodic%003i.pkl' % (path_to_filter, i), 'wb') as f:
        pickle.dump(argo_extract, f)


#  ----------------------------------------------------------------------------
def get_idx_from_tiles_lim(res, argodb):
    """Get the list of profile indices present in argodb that correspond
       to the list of wmos
       
       :rtype: dic"""
    if res['LONMIN_WITH_M'] > res['LONMAX_WITH_M']:
        idx = np.where((argodb['LATITUDE'] > res['LATMIN_WITH_M']) & (argodb['LATITUDE'] < res['LATMAX_WITH_M']) & ((argodb['LONGITUDE'] > res['LONMIN_WITH_M']) | (argodb['LONGITUDE'] < res['LONMAX_WITH_M'])))
    else:
        idx = np.where((argodb['LATITUDE'] > res['LATMIN_WITH_M']) & (argodb['LATITUDE'] < res['LATMAX_WITH_M']) & (argodb['LONGITUDE'] > res['LONMIN_WITH_M']) & (argodb['LONGITUDE'] < res['LONMAX_WITH_M']))
    argo_extract = extract_idx_from_argodb(argodb, idx)
    argo_extract['LATMIN_NO_M'] = res['LATMIN_NO_M']
    argo_extract['LATMAX_NO_M'] = res['LATMAX_NO_M']
    argo_extract['LONMIN_NO_M'] = res['LONMIN_NO_M']
    argo_extract['LONMAX_NO_M'] = res['LONMAX_NO_M']
    argo_extract['MARGINLAT'] = res['MARGINLAT']
    argo_extract['MARGINLON'] = res['MARGINLON']

    return argo_extract


#  ----------------------------------------------------------------------------
def get_idx_from_list_wmo(argodb, wmos):
    """Get the list of profile indices present in argodb that correspond
       to the list of wmos
       
       :rtype: list of int

    """
    infos = argotools.retrieve_infos_from_tag(argodb, argodb['TAG'])
    idx = []
    for w in wmos:
        idx += list(np.where(infos['WMO'] == w)[0])
    return idx


#  ----------------------------------------------------------------------------
def extract_idx_from_argodb(argodb, idx):
    """Return a argodb type dictionnary that is a subset of argodb and
       containing only entries given in idx (list)
       
       :rtype: dic

    """
    argodb_extract = {}
    for k in argodb.keys():
        argodb_extract[k] = argodb[k][idx]
    return argodb_extract


#  ----------------------------------------------------------------------------
def extract_idx_from_wmostats(wmostats, idx):
    """Return a wmostats type dictionnary that is a subset of wmostats and
       containing only entries given in idx (list)
       
       :rtype: dic

    """
    wmostats_extract = {}
    keys = wmostats.keys()
    keys.remove('N_WMO')
    for k in keys:
        wmostats_extract[k] = wmostats[k][idx]
    if type(idx) in [int, np.int64]:
        n_wmo = 1
    else:
        n_wmo = len(idx)
    wmostats_extract['N_WMO'] = n_wmo
    return wmostats_extract


def main():
    """Main function of stats.py"""
    creating_tiles()
    #  plot_map()

# ----------------------------------------
# Plot the tiles on map

def mbox(x1, x2, y1, y2, dlat, dlon, col, itile, m):
    m.plot([x1-dlon, x1-dlon, x2+dlon, x2+dlon, x1-dlon],
           [y1-dlat, y2+dlat, y2+dlat, y1-dlat, y1-dlat],
           color=col, latlon=True)
    x, y = m(0.5*(x1+x2), 0.5*(y1+y2-3))
    plt.text(x, y, '%i' % itile,
             color=col, horizontalalignment='center')


def plot_map():
    lat, lon, nlat, nlon, marginlat, marginlon = tile_definition()
    plt.figure(figsize=(12, 6))
    m = Basemap(projection='cea', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='coral', lake_color='white')

    itile = 0
    #  cols = plt.get_cmap('hot')
    for j in range(nlat):
        for i in range(nlon):
            itile = i + j*nlon
            if lon[i] - marginlon < -180:
                lon[i] = -178
            elif lon[i+1] + marginlon > 180:
                lon[i+1] = 178
            mbox(lon[i], lon[i+1], lat[j], lat[j+1], marginlat[j], marginlon,
                 'black', itile, m)

    plt.axis('tight')
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.savefig('tiles_summary.pdf')
    plt.show()


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    main()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
