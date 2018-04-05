# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:07:26 2018

@author: herry
"""

import research_tools as research
import argodb as argo
import argotools as argotools
import interpolation_tools as interpolation

zref = argotools.zref


#  ----------------------------------------------------------------------------
def creating_tiles(i):
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    wmodic = argo.read_wmodic()
    argodic = research.read_argo_filter(i)
    res = interpolation.interpolate_profiles(argodic, wmodic)
    print(res)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    creating_tiles(0)
