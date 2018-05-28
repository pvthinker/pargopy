#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:49:04 2018

@author: therry

Module contenant la série d'outil permettant la gestion d'un fichier netCDF :
    - Création du fichier, des variables, des dimensions, des attributs
    - Ecriture au sein du fichier de ces variables
    - Lecture du fichier

"""

#  from netCDF4 import Dataset
import os
import json as json
import datetime as datetime

#  ----------------------------------------------------------------------------
def create_dim(filename, zref, nlat, nlon, atlas_infos):
    """
    :param filename: Nom du fichier netCDF
    :param zref: Paliers de référence pour la profondeur
    :param nlon: Nombre de longitudes pour cette dalle
    :param nlat: Nombre de latitudes pour cette dalle
    :param atlas_infos: Dictionnaire contenant les champs suivants : (mode, date, reso, timeflag, typestat)
                      permettant de connaitre les informations sur les statistiques à calculer
                      ( on ne prendra ici que mode et date)
    
    Fonction utilisée pour créer les dimensions et les attributs d'un fichier 
    netCDF
    :rtype: None
    """

    creationdate = datetime.date.today()
    argodate = datetime.date(int(atlas_infos['DATE'][0]),
                             int(atlas_infos['DATE'][1]),
                             int(atlas_infos['DATE'][2]))

    print('create netCDF %s' % filename)

    with Dataset(filename, "w", format="NETCDF4") as rootgrp:

        rootgrp.createDimension('zref', len(zref))
        rootgrp.createDimension('lat', nlat)
        rootgrp.createDimension('lon', nlon)

        rootgrp.setncattr('using Argo profiles in data_mode', atlas_infos['MODE'])
        rootgrp.setncattr('using Argo profiles up to date', argodate.strftime('%Y-%m-%d'))
        rootgrp.setncattr('creation date', creationdate.strftime('%Y-%m-%d'))
        rootgrp.setncattr('author', 'Roullet, Guillaume')
        rootgrp.setncattr('email', 'roullet@univ-brest.fr')
        
#  ----------------------------------------------------------------------------
def create_var(filename, var_list):
    """
    :param filename: Nom du fichier netCDF
    :param var_list: Liste des variables à générer
    
    Create the netcdf file 'filename' for the list of variables 'var_list'
    the dimensions '(zref, lat, lon)' and the attributes of each variable is
    retrieved from the json file.
    
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
    :param filename: Nom du fichier netCDF
    :param var_list: Liste des variables à générer
    :param var_dic: Dictionnaire contenant les variables à écrire
    
    Write the list of variables 'var_list' in the netcdf file 'filename'. Variables are 
    transfered in the form of a dictionnary 'var_dic', where each entry is a numpy array.

    :rtype: None
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
    :param filename: Nom du fichier netCDF
    :param var_list: Liste des variables à générer
    
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

