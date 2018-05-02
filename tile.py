# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry

File used to generate the different tiles of the atlas

"""

import pickle
import time
import matplotlib.pyplot as plt
import argotools as argotools
import interpolation_tools as interpolation
import param as param
import numpy as np
import decorator as deco
import melted_functions as melted

zref = argotools.zref
path_to_tiles = param.path_to_tiles
path_to_data = param.path_to_data
path_to_filter = param.path_to_filter



#  ----------------------------------------------------------------------------
def creating_tiles(i):
    """Giving values to the variables
    
    :rtype: dic"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    #  wmodic = argotools.read_wmodic()
    #  argodic = argotools.read_argo_filter(i)
    wmodic = melted.read_dic('wmodic', path_to_data)
    argodic = melted.read_dic('argodic%003i' % i, path_to_filter)
    tile = interpolation.interpolate_profiles(argodic, wmodic)
    del(argodic, wmodic)
    tile['ZREF'] = zref

    #  write_tile(tile, i)
    melted.write_dic('tile%003i' % i, tile, path_to_tiles)

    return tile


#  ----------------------------------------------------------------------------
def write_tile(tile, i):
    """Write one tile in a .pkl file
    
    :rtype: None"""
    with open('%s/tile%003i.pkl' % (path_to_tiles, i), 'w') as f:
        pickle.dump(tile, f)


#  ----------------------------------------------------------------------------
def read_tile(i):
    """Read one of the tiles
    
    :rtype: dic"""
    print('read tile%003i.pkl' % i)
    with open('%s/tile%003i.pkl' % (path_to_tiles, i), 'r') as f:
        tile = pickle.load(f)
    return tile


#  ----------------------------------------------------------------------------
def plot_tile(i):
    """Plots the tiles values (Ti, Si, Ri) with the values non interpolate
    
    :rtype: None"""
    depth = zref
    #  argodic = argotools.read_argo_filter(i)
    #  tile = read_tile(i)
    argodic = melted.read_dic('argodic%003i' % i, path_to_filter)
    tile = melted.read_dic('tile%003i' % i, path_to_tiles)
    output = argotools.retrieve_infos_from_tag(argodic, argodic['TAG'])
    print(len(output['IDAC']))
    print(output['IDAC'])
    print(len(output['WMO']))
    print(output['WMO'])
    print(len(output['IPROF']))
    print(output['IPROF'])
    #  Can't find the path
    for i in range(len(output['WMO'])):
        dico = argotools.read_profile(output['IDAC'][0], output['WMO'][0], iprof=output['IPROF'][0], data=True)
    CT = tile['CT'][0]
    temp = dico['TEMP'][0]
    pres = dico['PRES'][0]
    plt.figure()
    plt.plot(CT, depth, '.')
    plt.plot(temp, pres, '+')
    plt.show()


def main(itile):
    """Main function of tile.py"""
    creating_tiles(itile)

#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
#==============================================================================
#     main(70)
#     main(71)
#==============================================================================
#      main(51)
#      main(50)
#      main(51)
    #  tiles = [131, 146, 151, 171, 190, 233, 234, 235, 244, 253, 254, 255, 264, 
    #         265, 272, 273, 274, 275, 276, 284, 285, 295, 296]
    #for t in tiles:
    main(189)
    deco.call_results('interpolation_tools.py')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
