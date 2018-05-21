#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:20:04 2018

@author: therry
"""
from netCDF4 import Dataset
import os
import json as json

typestat = ['general', 'zmean', 'zstd', 'zdz']

#  ----------------------------------------------------------------------------
def netCDF_dim_creation(filename, zref, nlat, nlon, mode, date):

    with Dataset(filename, "w", format="NETCDF4") as rootgrp:

        rootgrp.createDimension('zref', len(zref))
        rootgrp.createDimension('lat', nlat)
        rootgrp.createDimension('lon', nlon)
        rootgrp.setncattr('data_mode', mode)
        rootgrp.setncattr('year', date[0])
        rootgrp.setncattr('month', date[1])
        rootgrp.setncattr('day', date[2])


#  ----------------------------------------------------------------------------
def netCDF_var_creation(filename, var_list):
    """
    Create the netcdf file 'filename' for the list of variables 'var_list'
    the dimensions '(zref, lat, lon)' and the attributes of each variable is
    retrieved from the json file. The dimensions are automatically appended to
    the list of variables

    Warning 'filename' should be created first
    
    :rtype: None
    """

    var_list += ['zref', 'lon', 'lat']

    var_attributes = json.load(open('var_attributes.json'))

    if (os.path.isfile(filename)):
        print('%s exists' % filename)

        with Dataset(filename, "r+", format="NETCDF4") as rootgrp:

            for var in var_list:
                att = var_attributes[var]
                if att['dim'] == 1:
                    dimension = ('zref', )
                elif att['dim'] == 2:
                    dimension = ('lat', 'lon')
                elif att['dim']== 3:
                    dimension = ('zref', 'lat', 'lon')

                ncvar = rootgrp.createVariable(var, att['type'], dimension)
                ncvar.long_name = att['long_name']
                ncvar.units = att['unit']
    else:
        raise ValueError('%s should be created first' % filename)


#  ----------------------------------------------------------------------------
def netCDF_var_writing(filename, var_dic):
    """
    Write the list of variables in the netcdf file 'filename'. Variables are 
    transfered in the form of a dictionnary 'var_dic', where each entry is a numpy array.
    The dimensions '(zref, lat, lon)' should be in 'var_dic'.

    """

    if set(['zref', 'lon', 'lat']).issubset(set(var_dic.keys())):
        pass
    else:
        raise ValueError('zref, lon and lat should be in var_dic')

    if (os.path.isfile(filename)):
    
        with Dataset(filename, "r+", format="NETCDF4") as f:

            for var in var_dic.keys():
                var_data = var_dic[var]
                ndims = len(var_data.shape)
                if ndims == 1:
                    f.variables[var][:] = var_data
                elif ndims == 2:
                    f.variables[var][:, :] = var_data
                elif ndims == 3:
                    f.variables[var][:, :, :] = var_data

    else:
        raise ValueError('%s should be created first' % filename)

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

