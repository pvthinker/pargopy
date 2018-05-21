#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:20:04 2018

@author: therry
"""
from netCDF4 import Dataset
import os
import json as json
import datetime as datetime

typestat = ['general', 'zmean', 'zstd', 'zdz']

#  ----------------------------------------------------------------------------
def create_dim(filename, zref, nlat, nlon, mode, date):

    creationdate = datetime.date.today()
    argodate = datetime.date(int(date[0]),int(date[1]),int(date[2]))

    print('create netCDF %s' % filename)

    with Dataset(filename, "w", format="NETCDF4") as rootgrp:

        rootgrp.createDimension('zref', len(zref))
        rootgrp.createDimension('lat', nlat)
        rootgrp.createDimension('lon', nlon)

        rootgrp.setncattr('using Argo profiles in data_mode', mode)
        rootgrp.setncattr('using Argo profiles up to date', argodate.strftime('%Y-%m-%d'))
        rootgrp.setncattr('creation date', creationdate.strftime('%Y-%m-%d'))
        rootgrp.setncattr('author', 'Roullet, Guillaume')
        rootgrp.setncattr('email', 'roullet@univ-brest.fr')
        
#  ----------------------------------------------------------------------------
def create_var(filename, var_list):
    """
    Create the netcdf file 'filename' for the list of variables 'var_list'
    the dimensions '(zref, lat, lon)' and the attributes of each variable is
    retrieved from the json file.

    Warning 'filename' should be created first
    
    :rtype: None
    """

    var_attributes = json.load(open('var_attributes.json'))

    if (os.path.isfile(filename)):

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
def write_var(filename, var_list, var_dic):
    """
    Write the list of variables 'var_list' in the netcdf file 'filename'. Variables are 
    transfered in the form of a dictionnary 'var_dic', where each entry is a numpy array.

    """

    if any([var_dic.has_key(v) for v in var_list]):
        pass
    else:
        raise ValueError('some elements in var_list are missing in var_dic')

    if (os.path.isfile(filename)):
    
        with Dataset(filename, "r+", format="NETCDF4") as f:

            for var in var_list:
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
def read_var(filename, var_list):
    """
    Read the list of variables 'var_list' from netcdf file 'filename'
    Return the result as a dictionnary of numpy arrays

    :rtype: dict
    """

    var_dic = {}
    if (os.path.isfile(filename)):        

        var_attributes = json.load(open('var_attributes.json'))

        with Dataset(filename, "r", format="NETCDF4") as f:
            for var in var_list:
                att = var_attributes[var]
                if att['dim'] == 1:
                    var_dic[var] = f.variables[var][:]
                elif att['dim'] == 2:
                    var_dic[var] = f.variables[var][:, :]
                elif att['dim']== 3:
                    var_dic[var] = f.variables[var][:, :, :]

    else:
        print('Warning: %s does not exist' % filename)

    return var_dic

