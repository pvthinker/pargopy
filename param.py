# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:01:25 2018

@author: herry
"""
#  location is used to know where you are working
#  if you're working on datarmor, location takes the value 'DATARMOR'
#  if you're working on herry-s local, location takes the value 'HERRY'
location = 'HERRY'
#  location = 'DATARMOR'

if location == 'DATARMOR':
    path_to_argo = '/home/oo26/coriolis/co05/co0508/dac'

    path_to_pargopy = 'home2/datahome/therry/pargopy/'

    path_to_data = 'home2/datawork/therry/data'

    path_to_filter = 'home2/datawork/therry/filter'

    path_to_stats = 'home2/datawork/therry/stats'

    path_to_tiles = 'home2/datawork/therry/tiles'

elif location == 'HERRY':
    path_to_argo = '/net/alpha/exports/sciences/data/ARGO/ARGO/201602-ArgoData'

    path_to_pargopy = '/home/herry/Documents/pargopy/'

    path_to_data = '/local/tmp/1/herry/pargopy/data'

    path_to_filter = '/local/tmp/1/herry/pargopy/filter'

    path_to_stats = '/local/tmp/1/herry/pargopy/stats'

    path_to_tiles = '/local/tmp/1/herry/pargopy/tiles'
