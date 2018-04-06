# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry
"""

import pickle
import time
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

    write_tile(tile, i)

    return tile


#  ----------------------------------------------------------------------------
def write_tile(tile, i):
    with open('%s/tile%003i.pkl' % (path_localdata, i), 'w') as f:
        pickle.dump(tile, f)


#  ----------------------------------------------------------------------------
def read_argo_filter(i):
    print('read tile%003i.pkl' % i)
    with open('%s/tile%003i.pkl' % (path_localdata, i), 'r') as f:
        tile = pickle.load(f)
    return tile


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    creating_tiles(0)
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
