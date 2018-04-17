# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry
"""

import pickle
import time
import matplotlib.pyplot as plt
import research_tools as research
import argotools as argotools
import interpolation_tools as interpolation
import param as param

zref = argotools.zref
path_localdata = param.path_to_tiles


#  ----------------------------------------------------------------------------
def creating_tiles(i):
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    wmodic = argotools.read_wmodic()
    argodic = argotools.read_argo_filter(i)
    print('argodic red')
    tile = interpolation.interpolate_profiles(argodic, wmodic)
    tile['ZREF'] = zref

    write_tile(tile, i)

    return tile


#  ----------------------------------------------------------------------------
def write_tile(tile, i):
    with open('%s/tile%003i.pkl' % (path_localdata, i), 'w') as f:
        pickle.dump(tile, f)


#  ----------------------------------------------------------------------------
def read_tile(i):
    print('read tile%003i.pkl' % i)
    with open('%s/tile%003i.pkl' % (path_localdata, i), 'r') as f:
        tile = pickle.load(f)
    return tile


#  ----------------------------------------------------------------------------
def plot_tile(i):
    """Plots the tiles values (Ti, Si, Ri) with the values non interpolate"""
    depth = zref
    argodic = argotools.read_argo_filter(i)
    tile = read_tile(i)
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
    print('Main of tile')
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
    main(299)
    #  73, 41, 118
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
