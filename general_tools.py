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

import jdcal
import numpy as np

import param as param


def unmask(data):
    """ transform masked array into regular numpy array """
    data_out = {}
    for k in data.keys():
        if type(data[k]) is np.ma.core.MaskedArray:
            data_out[k] = data[k].data
        else:
            data_out[k] = data[k]
    return data_out


def bytes2str(data):
    """ byte strings into strings"""
    data_out = {}
    for k in data.keys():
        data_out[k] = data[k]
        if type(data[k]) is np.ndarray:
            firstelem = data_out[k].ravel()[0]
            if type(firstelem) is np.bytes_:
                data_out[k] = np.asarray(data[k].data, dtype=str)
    return data_out


def get_tag(kdac, wmo, kprof):
    """Compute the tag number of a profile

    The inverse of get_tag() is retrieve_infos_from_tag()

    :rtype: int
    """
    if type(kprof) == int:
        kprof = [kprof]

    elif hasattr(kprof, '__iter__'):
        pass

    else:
        raise ValueError('kprof must an int, a float or be iterable')

    if max(kprof) > 1000:
        raise ValueError("kprof > 1000, tags are no longer unique")

    tag = [(kdac*10000000+wmo)*1000+k for k in kprof]

    if len(kprof) == 1:
        tag = tag[0]

    return tag


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


def dac_from_wmo(wmodic, wmo):
    """
    Retrieve the dac of a wmo
    :rtype: list of dac
    """

    dac = ''
    for d in param.daclist:
        if wmo in wmodic[d]:
            dac = d
    return dac


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


def count_profiles_in_database(wmostats):
    """
    Count the total number of profiles in database
    :rtype: int
    """

    nbprofiles = 0
    for nbpr in wmostats['N_PROF']:
        nbprofiles += nbpr
    print('number of new profiles in Argo database: %i' % nbprofiles)
    return nbprofiles


def compute_weight(x, y, lon, lat, reso):
    """Compute the weight between points (x, y) and point (lon, lat) with
    a gaussian filter """
    dist = dist_sphe(x, y, lon, lat)
    weight = np.exp(-0.5*(dist/reso)**2)
    return weight


def dist_sphe(x, y, lon, lat):
    """Compute the spherical arc between two points on the unit sphere"""
    return np.arccos(np.sin(lat)*np.sin(y)+np.cos(lat)*np.cos(y)*np.cos(lon-x))


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
