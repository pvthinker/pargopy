#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:50:24 2018

@author: therry
"""

import tile as ti
import param as param
import general_tools as tools
import database as db

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend("TkAgg")

def plot_figures():
    """
    Fonction utilisée pour afficher les et enregistrer les figures
    """
    pass


itile = 32

tile = ti.read_argo_tile(itile)
zref_profile = ti.read_zref_profiles(itile)
var = zref_profile['CT']

lat = tile['LATITUDE'].loc[var.index]
lon = tile['LONGITUDE'].loc[var.index]

#==============================================================================
# plt.figure()
# plt.title('Temperature interpolee a %1f metres de profondeur' % var.columns[0])
# plt.xlabel('LONGITUDE')
# plt.ylabel('LATITUDE')
# plt.scatter(lon, lat, c=var.iloc[:,0], s=10, cmap='jet')
# plt.colorbar()
# plt.show()
#==============================================================================

wmos = []

tag = var.index[0]
tag_infos = tools.retrieve_infos_from_tag(tag)

res = db.read_profile(param.daclist[tag_infos['IDAC']], tag_infos['WMO'], iprof=tag_infos['IPROF'], data=True)

plt.figure(1)
plt.title('Temperature d un profil apres interpolation')
plt.xlabel('Temperature')
plt.ylabel('Profondeur')
plt.scatter(res['TEMP'], -res['PRES'], s=10, cmap='red', label='Temperature mesuree')
plt.scatter(var.iloc[0], -param.zref, s=10, label='Temperature interpolee')
plt.legend(loc=2)

#  plt.show()

#==============================================================================
# plt.figure(2)
# plt.title('Temperature d un profil avant interpolation')
# plt.xlabel('Temperature')
# plt.ylabel('Profondeur')
# plt.scatter(res['TEMP'], -res['PRES'], s=10)
#==============================================================================
plt.show()