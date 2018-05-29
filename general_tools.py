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
    return juliandayf - 2433282.5