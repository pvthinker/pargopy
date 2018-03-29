# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:10:24 2018

@author: herry
"""
import jdcal
import os
import numpy as np
import param as param
from netCDF4 import Dataset


#  ----------------------------------------------------------------------------
def path_to_wmo():
    """Method creating the first folder containing all the informations of
    the differents profiles known on Argo"""
    path_to_wmo, path_to_argolog = param.data_argotools()
    #  Creating a list with the name of the following repository
    list_folder = []
    list_folder += os.listdir(path_to_wmo)

    #  Sorting the results of list_folder_0 to keep only the good results (without .tar.gz, ...)
    true_dir = [i for i in list_folder if i[-8:] == 'ArgoData']

    for i in range(len(true_dir)):
        list_wmoid_bis = [[], [], [], [], [], [], [], [], [], [], []]
        len_argolog = 2000000
        list_wmoid = np.zeros((len_argolog,), dtype=int)
        size = 0
        list_folder_1 = []
        first_dir_list = []
        dir_file_list = []
        final_dir_list = []
        dac_list = []
        wmo_list = []
        list_folder_1 += os.listdir('%s/%s' % (path_to_wmo, true_dir[i]))
        #  Listing the following repository with the entire adress found
        #  with the first_dir
        current_directory = '%s/%s' % (path_to_wmo, true_dir[i])
        first_dir_list.append(current_directory)

        #  Listing the following datacenter (aoml, bodc, coriolis, ...)

        dac_list += [i for i in list_folder_1 if len(i) <= 8 and i[:] != 'doc']

    for i in range(len(dac_list)):
        current_directory_2 = '%s/%s' % (current_directory, dac_list[i])
        #  list_wmoid[size:size+len(current_directory_2)] = os.listdir('%s' % (current_directory_2))
        #  size += len(current_directory_2)
        list_wmoid_bis[i] += os.listdir('%s' % (current_directory_2))

        for j in range(len(list_wmoid_bis[i])):
            current_directory_3 = "%s/%s" % (current_directory_2, list_wmoid_bis[i][j])
            final_dir_list.append(current_directory_3)
            directory_and_file = "%s/%s/%s" % (current_directory_2, list_wmoid_bis[i][j], list_wmoid_bis[i][j])
            dir_file_list.append(directory_and_file)

    for i, wmo in enumerate(list_wmoid_bis):
        list_wmoid[size:size+len(wmo)] = wmo
        size += len(wmo)

    return (dac_list, wmo_list, final_dir_list, list_wmoid, dir_file_list)

#  ----------------------------------------------------------------------------
    def reading_variables(minlat, maxlat, minlon, maxlon):
        """Salvation of the profiles datas"""
        #  Reading the netCDF file fom ARGO
        if (os.path.isfile(path_to_argolog)):
            f = Dataset(path_to_argolog, "r", format="NETCDF4")
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            idx = np.where((lats > minlat) & (lats < maxlat) & (lons > minlon) & (lons < maxlon))
            print(len(idx))
            latitudes = f.variables['lat'][idx]
            longitudes = f.variables['lon'][idx]
            prof = f.variables['prof'][idx]
            update = f.variables['update'][idx]
            juld = f.variables['juld'][idx]
            dac = f.variables['dac'][idx]
            wmo = f.variables['wmo'][idx]

            f.close()

        return (latitudes, longitudes, prof, update, juld, dac, wmo)

#  ----------------------------------------------------------------------------
def conversion_juld_gregd(juld):
    """Method converting julian day into gregorian day"""
    #  lats, lons, juld = self.reading_variables()
    #  2433282.5000000 corresponds to the Argo origin date
    gregday = jdcal.jd2jcal(2433282.500000, juld)
    print('This Julian Day corresponds to {0}/{1}/{2}'.format(gregday[2], gregday[1], gregday[0]))
    return(gregday)


#  ----------------------------------------------------------------------------
def conversion_gregd_juld(day, month, year):
    """Method converting gregorian day into julian day"""
    #  Petit décalage possible
    julianday = jdcal.gcal2jd(year, month, day)
    print('This Gregorian Day corresponds to {0}'.format((julianday[0]+julianday[1]+0.5)-2433282.500000))
