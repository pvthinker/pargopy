# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry

File used to generate the different tiles of the atlas

"""

import time
import matplotlib.pyplot as plt
import argotools as argotools
import interpolation_tools as interpolation
import param as param
#import decorator as deco

zref = argotools.zref
path_to_tiles = param.path_to_tiles
path_to_data = param.path_to_data
path_to_filter = param.path_to_filter


#  ----------------------------------------------------------------------------
def generate_tile(i):
    """Interpolate all Argo profiles in tile 'i' onto 'zref' depths. Save
    the result in the 'tile%003i.pkl' file

    """

    wmodb = argotools.read_dic('wmodic', path_to_data)
    argodb = argotools.read_dic('argo%003i' % i, path_to_filter)
    zrefprofiles = interpolation.interpolate_profiles(argodb, wmodb)
    zrefprofiles['ZREF'] = zref

    argotools.write_dic('tile%003i' % i, zrefprofiles, path_to_tiles)

    # try to reduce memory leakage when processessing all the tiles
    del(argodb, wmodb)


#  ----------------------------------------------------------------------------
def generate_argotiles():
    """Generate the argodb dictionnary for each tile and save it in
    'argo%003i'
    """

    argodb = argotools.read_dic('argodb', path_to_data)
    lat, lon, nlat, nlon, marginlat, marginlon = argotools.tile_definition()
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

            argo_extract = argotools.extract_idx_inside_tile(res, argodb)
            argotools.test_tiles(argo_extract, k)
            argotools.write_dic('argo%003i' %
                                k, argo_extract, path_to_filter)
            k += 1


#  ----------------------------------------------------------------------------
def plot_tile(i):
    """Plots the tiles values (Ti, Si, Ri) with the values non interpolate

    :rtype: None"""
    depth = zref
    argo = argotools.read_dic('argo%003i' % i, path_to_filter)
    tile = argotools.read_dic('tile%003i' % i, path_to_tiles)
    output = argotools.retrieve_infos_from_tag(argo['TAG'])
    print(len(output['IDAC']))
    print(output['IDAC'])
    print(len(output['WMO']))
    print(output['WMO'])
    print(len(output['IPROF']))
    print(output['IPROF'])
    #  Can't find the path
    for i in range(len(output['WMO'])):
        dico = argotools.read_profile(
            output['IDAC'][0], output['WMO'][0],
            iprof=output['IPROF'][0], data=True)
    CT = tile['CT'][0]
    temp = dico['TEMP'][0]
    pres = dico['PRES'][0]
    plt.figure()
    plt.plot(CT, depth, '.')
    plt.plot(temp, pres, '+')
    plt.show()


def main(itile):
    """Main function of tile.py"""
    generate_tile(itile)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
# ==============================================================================
#     main(70)
#     main(71)
# ==============================================================================
#      main(51)
#      main(50)
#      main(51)
    #  tiles = [131, 146, 151, 171, 190, 233, 234, 235, 244, 253, 254, 255, 264,
    #         265, 272, 273, 274, 275, 276, 284, 285, 295, 296]
    # for t in tiles:
    main(189)
    deco.call_results('interpolation_tools.py')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
