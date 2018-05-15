#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:20:04 2018

@author: therry
"""
from netCDF4 import Dataset
import os
import simplejson as json

typestat = ['general', 'zmean', 'zstd', 'zdz']

#  ----------------------------------------------------------------------------
def netCDF_dim_creation(filename, zref, nlat, nlon, mode, date):

    rootgrp = Dataset(filename, "w", format="NETCDF4")

    rootgrp.createDimension('zref', len(zref))
    rootgrp.createDimension('lat', nlat)
    rootgrp.createDimension('lon', nlon)
    rootgrp.setncattr('data_mode', mode)
    rootgrp.setncattr('year', date[0])
    rootgrp.setncattr('month', date[1])
    rootgrp.setncattr('day', date[2])

    rootgrp.close()


#  ----------------------------------------------------------------------------
def netCDF_var_creation(filename, var_choice):
    """
    Function creating the netDCF variables.
    This function takes its values from a .json file where all the desired
    values are handly referenced
    WARNING : You must open and close the ncfile before and after calling
              this function, the ame is to give the name you want to the ncfile
    
    :rtype: None
    """

    if (os.path.isfile(filename)):
        print('filename exists')
        rootgrp = Dataset(filename, "r+", format="NETCDF4")

        data = json.load(open('pargopy_var.json'))
        nc_var = {}

        for g in data['general']:
            if g['name'] in var_choice:
                pass
            else:
                var_choice.append(g['name'])
            print(var_choice)

        for v in var_choice:
            for t in typestat:
                for d in data[t]:
                    if v == d['name']:
                        if d['dim'] == 1:
                            dimension = ('zref', )
                        elif d['dim'] == 2:
                            dimension = ('lat', 'lon')
                        elif d['dim']== 3:
                            dimension = ('zref', 'lat', 'lon')

                        if d['name'] in nc_var:
                            pass
                        else:
                            nc_var[d['name']] = rootgrp.createVariable(d['name'], d['type'], dimension)
                            nc_var[d['name']].long_name = d['long_name']
                            nc_var[d['name']].units = d['unit']

        rootgrp.close()

#  ----------------------------------------------------------------------------
def netCDF_var_writing(filename, var_choice, res):
    """
    Function writing the netDCF variables.
    This function takes its values from a .json file where all the desired
    values are handly referenced and from res(type = dict)
    WARNING : You must open and close the ncfile before and after calling
              this function, the ame is to give the name you want to the ncfile
    
    :rtype: None
    """

    if (os.path.isfile(filename)):
        print('filename exists')
        print(filename)
        print(var_choice)
        f = Dataset(filename, "r+", format="NETCDF4")

        data = json.load(open('pargopy_var.json'))

        for g in data['general']:
            if g['name'] in var_choice:
                pass
            else:
                var_choice.append(g['name'])
            print(var_choice)

        for v in var_choice:
            for t in typestat:
                for d in data[t]:
                    if v == d['name']:
                        #  print(d['name'])
                        if d['dim'] == 1:
                            f.variables[d['name']][:] = res[d['name']]
                        elif d['dim'] == 2:
                            f.variables[d['name']][:, :] = res[d['name']]
                        elif d['dim']== 3:
                            f.variables[d['name']][:, :, :] = res[d['name']]

        f.close()


#  ----------------------------------------------------------------------------
def netCDF_var_reading(filename, var_choice):
    """
    Function reading the netDCF files.
    This function takes its values from a .json file where all the desired
    values are handly referenced and from res(type = dict)
    WARNING : You must open and close the ncfile before and after calling
              this function, the ame is to give the name you want to the ncfile
    
    :rtype: dict
    """

    if (os.path.isfile(filename)):
        f = Dataset(filename, "r", format="NETCDF4")

        data = json.load(open('pargopy_var.json'))
        res = {}

        for g in data['general']:
            var_choice.append(g['name'])

        for v in var_choice:
            for t in typestat:
                for d in data[t]:
                    if v == d['name']:
                        if d['dim'] == 1:
                            res[d['name']] = f.variables[d['name']][:]
                        elif d['dim'] == 2:
                            res[d['name']] = f.variables[d['name']][:, :]
                        elif d['dim']== 3:
                            res[d['name']] = f.variables[d['name']][:, :, :]

        f.close()
    else:
        res = {}

    return res

