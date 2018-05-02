#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:20:04 2018

@author: therry
"""
from netCDF4 import Dataset
#  import json as json
import simplejson as json
import argotools as argotools

zref = argotools.zref

def netCDF_var_creation(typestat):
    """
    Function creating the netDCF variables.
    This function takes its values from a .json file where all the desired
    values are handly referenced
    """
    data = json.load(open('%s_var.json' % typestat))

    rootgrp = Dataset('test.nc', "w", format="NETCDF4")
    nlat = 15
    nlon = 20
    rootgrp.createDimension('zref', len(zref))
    rootgrp.createDimension('lat', nlat)
    rootgrp.createDimension('lon', nlon)
    dimension = ('lat', 'lon')
    for d in data[typestat]:
        print(d)
        if d['dim'] == 1:
            dimension = ('zref', )
        elif d['dim'] == 2:
            dimension = ('lat', 'lon')
        elif d['dim']== 3:
            dimension = ('zref', 'lat', 'lon')
        
        data[d['name']] = rootgrp.createVariable(d['name'], d['type'], dimension)
        data[d['name']].long_name = d['long_name']
        data[d['name']].units = d['unit']

    rootgrp.close()


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    netCDF_var_creation('zstd')

