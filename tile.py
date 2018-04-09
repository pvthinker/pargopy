# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry
"""

import pickle
import time
import matplotlib.pyplot as plt
import research_tools as research
import argodb as argo
import argotools as argotools
import interpolation_tools as interpolation

zref = argotools.zref
path_localdata = '/local/tmp/1/herry/pargopy/tiles/'


#  ----------------------------------------------------------------------------
def creating_tiles(i):
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    wmodic = argo.read_wmodic()
    argodic = research.read_argo_filter(i)
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
    argodic = research.read_argo_filter(i)
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



#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    creating_tiles(0)
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
