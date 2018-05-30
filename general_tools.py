#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:50:19 2018

@author: therry

Module contenant les outils utilisés dans la plupart des modules:
    - Lecture/Ecriture au sein de fichiers pickle
    - Outil de navigation au sein d'une DataFrame
    - ...

"""

import os
import jdcal
from netCDF4 import Dataset
import numpy as np

import param as param


def read_profile(dac, wmo, iprof=None,
                 header=False, data=False,
                 headerqc=False, dataqc=False,
                 verbose=True):
    """
    :param dac: DAC du profil recherché
    :param wmo: WMO du profil recherché
    :param iprof: Numéro du profil recherché
    :param header: Sélectionne seulement LATITUDE, LONGITUDE et JULD
    :param headerqc: Sélectionne seulement POSITION_QC et JULD_QC
    :param data: Sélectionne TEMP, PSAL et PRES
    :param dataqc: Sélectionne TEMP_QC, PSAL_QC et PRES_QC
    :param verbose: ???
    
    Les valeurs sélectionnée grâce aux arguments passés à la fonction définissent
    la DataFrame que retournera celle-ci.
    
    Basic driver to read the \*_prof.nc data file

    The output is a dictionnary of vectors
    - read one or all profiles read the header (lat, lon, juld) or not
    - read the data or not always return IDAC, WMO, N_PROF, N_LEVELS
    - and DATA_UPDATE (all 5 are int)

    :rtype: dict
    """

    key_header = ['LATITUDE', 'LONGITUDE', 'JULD']
    key_headerqc = ['POSITION_QC', 'JULD_QC']
    key_data = ['TEMP', 'PSAL', 'PRES']
    key_dataqc = ['TEMP_QC', 'PSAL_QC', 'PRES_QC']

    if type(dac) is int:
        dac = param.daclist[dac]

    filename = param.get_argo_filename(dac, wmo)

    if verbose:
        print('/'.join(filename.split('/')[-3:]))

    argo_profile_dic = {}

    required_keys = set(['TEMP', 'PSAL', 'PRES'])

    if (os.path.isfile(filename)):
        with Dataset(filename, "r", format="NETCDF4") as f:
            argo_profile_dic['DACID'] = param.daclist.index(dac)
            argo_profile_dic['WMO'] = wmo
            argo_profile_dic['N_PROF'] = len(f.dimensions['N_PROF'])
            argo_profile_dic['N_LEVELS'] = len(f.dimensions['N_LEVELS'])
            # DATE_UPDATE is an array of 14 characters in the *_prof.nc
            # we transform it into an int
            # YYYYMMDDhhmmss
            argo_profile_dic['DATE_UPDATE'] = ''.join(f.variables['DATE_UPDATE'][:])

            keyvar = set(f.variables.keys())

            if required_keys.issubset(keyvar):
                argo_profile_dic['TSP_QC'] = '1'
            else:
                argo_profile_dic['TSP_QC'] = '2'

            if header or headerqc or data or dataqc:
                if iprof is None:
                    idx = range(argo_profile_dic['N_PROF'])
                    argo_profile_dic['IPROF'] = np.arange(argo_profile_dic['N_PROF'])
                else:
                    idx = iprof
                    argo_profile_dic['IPROF'] = iprof

            if header:
                for key in key_header:
                    argo_profile_dic[key] = f.variables[key][idx]
                    argo_profile_dic['DATA_MODE'] = np.asarray(
                        [c for c in f.variables['DATA_MODE'][idx]])

            if headerqc:
                for key in key_headerqc:
                    argo_profile_dic[key] = f.variables[key][idx]

            if data:
                for key in key_data:
                    if argo_profile_dic['TSP_QC'] == '1':
                        argo_profile_dic[key] = f.variables[key][idx, :]
                    else:
                        argo_profile_dic[key] = np.NaN+np.zeros(
                            (argo_profile_dic['N_PROF'], argo_profile_dic['N_LEVELS']))

            if dataqc:
                for key in key_dataqc:
                    if argo_profile_dic['TSP_QC'] == '1':
                        argo_profile_dic[key] = f.variables[key][idx]
                    else:
                        argo_profile_dic[key] = np.zeros(
                            (argo_profile_dic['N_PROF'], argo_profile_dic['N_LEVELS']), dtype=str)

    return argo_profile_dic


def get_tag(kdac, wmo, kprof):
    """
    :param kdac: Index of the dac (aoml = 1, bodc = 2, coriolis = 3, ...)
    :param wmo: WMO number
    :param kprof: Index of the profile
    
    Compute the tag number of a profile

    The inverse of get_tag() is retrieve_infos_from_tag()

    :rtype: int
    """
    if kprof > 1000:
        raise ValueError("kprof > 1000, the tag may be wrong")

    return (kdac*10000000+wmo)*1000+kprof


def retrieve_infos_from_tag(tag):
    """
    :param tag: tag on which you are looking for informations

    Retrieve idac, wmo and iprof from tag (array of int)

    It is the inverse of get_tag()

    :rtype: dic
    """

    iprof = tag % 1000
    tag = (tag-iprof) // 1000
    wmo = tag % 10000000
    tag = (tag-wmo) // 10000000
    idac = tag
    tag_infos = {'IDAC': idac, 'WMO': wmo, 'IPROF': iprof}
    return tag_infos


def ij2tile(i, j):
    """
    :param i: Longitude du point étudié
    :param j: Latitude du point étudié
    
    Fonction retournant le numéro de la dalle à laquelle appartient un point
    dont les coordonnées sont i et j (longitude, latitude)
    TO DO: Changer le 20 par une valeur (nlat) pour rendre le code plus flexible
    
    :rtype: int
    """

    return i + j*20


#  ----------------------------------------------------------------------------
def conversion_juld_gregd(juld):
    """
    :param juld: Date en calendrier julien à convertir
    
    Fonction convertissant un julian day en gregorian day

    :rtype: list of int"""

    gregday = jdcal.jd2gcal(2433282.5, juld)
    print('This Julian Day corresponds to {0}/{1}/{2}'.format(
        gregday[2], gregday[1], gregday[0]))

    return(gregday)


#  ----------------------------------------------------------------------------
def conversion_gregd_juld(date):
    """
    :param date: Liste (année, mois, jour) contenant une date en calendrier
                 grégorien
    
    Fonction convertissant une date du calendrier grégorien en un julian day

    :rtype: float"""

    julianday = jdcal.gcal2jd(date[0], date[1], date[2])
    juliandayf = julianday[0] + julianday[1]
    return juliandayf - 2433282