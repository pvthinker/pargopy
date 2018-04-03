# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:22:03 2018

@author: herry
"""

from netCDF4 import Dataset
import numpy as np
import os
import math
import time
import numpy as np
import argotools as argotools
import param as param
tmps1 = time.time()


class ArgoDB:
    """Class to test netCDF files creation"""

    def __init__(self):
        """Initialization of the table"""
        #  Creating and opening the netCDF file
        self.path_to_argolog = param.data_argodb
        self.rootgrp = Dataset(self.path_to_argolog, "w", format="NETCDF4")
        print self.rootgrp.data_model
        self.variables_creation()
        self.writing_variables()

#  ----------------------------------------------------------------------------
    def variables_creation(self):
        """Creating the variables"""
        #  Generation of the dimension of the table

        line_number = self.rootgrp.createDimension('line_number', None)

        #  Generation of the variables of the table

        self.update = self.rootgrp.createVariable("update", "d", ('line_number', ))
        self.update.long_name = 'last update date'
        self.update.units = 'datetime'

        self.juliandays = self.rootgrp.createVariable("juld", "f4", ('line_number', ))
        self.juliandays.long_name = 'date'
        self.juliandays.units = 'julian day'

        self.prof = self.rootgrp.createVariable('prof', 'i2', ('line_number', ))
        self.prof.long_name = 'index of the profile in the profiles dir'
        self.prof.units = 'dimensionless'

        self.latitudes = self.rootgrp.createVariable("lat", "f4", ('line_number', ))
        self.latitudes.long_name = 'latitude'
        self.latitudes.units = 'degrees'

        self.longitudes = self.rootgrp.createVariable("lon", "f4", ('line_number', ))
        self.longitudes.long_name = 'longitude'
        self.longitudes.units = 'degrees'

        self.dac = self.rootgrp.createVariable('dac', 'b', ('line_number', ))
        self.dac.long_name = 'index of the dac'
        self.dac.units = 'dimensionless'

        self.wmo = self.rootgrp.createVariable('wmo', 'i', ('line_number', ))
        self.wmo.long_name = 'wmoid of the wmo'
        self.wmo.units = 'dimensionless'

        self.flag = self.rootgrp.createVariable('flag', 'i', ('line_number', ))
        self.flag.long_name = 'flag to know if the profile is OK'
        self.flag.units = 'dimensionless'

        return (self.latitudes, self.longitudes, self.juliandays, self.dac,
                self.wmo, self.update, self.prof, self.flag)

#  ----------------------------------------------------------------------------
    def write_profiles_notice(self, flag):
        """Method creating the .txt file with the following profiles
        informations"""

        d_list, dir_list, matrix_wmoid, list_wmoid, file_name = path_to_wmo()

        for i in range(len(dir_list)):
            fichier = open("/local/tmp/1/herry/profiles_notice.txt", "w")
            fichier.write('%s\n' % (flag[i]))
            fichier.close()

#  ----------------------------------------------------------------------------
    def reading_variables(self):
        """Salvation of the profiles datas"""
        #  Reading the netCDF file fom ARGO
        d_list, dir_list, matrix_wmoid, list_wmoid, file_name_list = path_to_wmo()
        len_argolog = 2000000
        lats1 = np.zeros((len_argolog,))
        lons1 = np.zeros((len_argolog,))
        juld1 = np.zeros((len_argolog,))
        dac1 = np.zeros((len_argolog,), dtype = int)
        wmo1 = np.zeros((len_argolog,), dtype = int)
        update1 = np.zeros((len_argolog,), dtype = int)
        prof = np.zeros((len_argolog,), dtype = int)
        flag1 = np.zeros((len_argolog,),dtype = int)
        size = 0

        for i, fl in enumerate(file_name_list):
            if (os.path.isfile('%s_prof.nc' % (fl))):
                f = Dataset('%s_prof.nc' % (fl), "r", format="NETCDF4")
                lats = f.variables['LATITUDE'][0:]
                lons = f.variables['LONGITUDE'][0:]
                juld = f.variables['JULD'][0:]
                update = f.variables['DATE_UPDATE'][:]
                lats1[size:size+len(lats)] = lats
                lons1[size:size+len(lats)] = lons
                juld1[size:size+len(lats)] = juld
                dac1 = self.finding_dac(fl, dac1, lats, size)
                wmo1 = self.finding_wmo(fl, wmo1, list_wmoid, lats, i, size)
                update1 = self.finding_update_date(fl, update, update1, lats, size)
                prof = self.finding_prof(prof, lats, size)
                flag = self.creating_flags(f, flag1, lats, size)
                size += len(lats)
                print('Apres passage dans finding_wmo numéro %i' % i)

                f.close()
        print('Quantity of values found : %i' % size)

        return lats1, lons1, juld1, dac1, wmo1, update1, prof, flag, size

#  ----------------------------------------------------------------------------
    def finding_dac(self, fl, dac1, lats, size):
        """Giving values to the dac"""
        #  Matching the dac with the following values of the netCDF doc
        if fl.find('aoml') != -1:
            dac1[size:size+len(lats)] = 0
        elif fl.find('bodc') != -1:
            dac1[size:size+len(lats)] = 1
        elif fl.find('coriolis') != -1:
            dac1[size:size+len(lats)] = 2
        elif fl.find('csio') != -1:
            dac1[size:size+len(lats)] = 3
        elif fl.find('csiro') != -1:
            dac1[size:size+len(lats)] = 4
        elif fl.find('incois') != -1:
            dac1[size:size+len(lats)] = 5
        elif fl.find('jma') != -1:
            dac1[size:size+len(lats)] = 6
        elif fl.find('kma') != -1:
            dac1[size:size+len(lats)] = 7
        elif fl.find('kordi') != -1:
            dac1[size:size+len(lats)] = 8
        elif fl.find('meds') != -1:
            dac1[size:size+len(lats)] = 9
        elif fl.find('nmdis') != -1:
            dac1[size:size+len(lats)] = 10

        return(dac1)

#  ----------------------------------------------------------------------------
    def finding_wmo(self, fl, wmo1, list_wmoid, lats, i, size):
        """Giving values to the wmo"""
        #  Matching the wmo with the following values of the netCDF doc
        if fl.find(str(list_wmoid[i])) != -1:
            wmo1[size:size+len(lats)] = list_wmoid[i]
            print(list_wmoid[i])

        return(wmo1)

#  ----------------------------------------------------------------------------
    def finding_prof(self, prof, lats, size):
        """Giving values to the dac"""
        #  Matching the profile rank with the following values of the netCDF doc
        i = 0
        while i < len(lats):
            prof[size+i] = i
            i = i + 1

        return(prof)

#  ----------------------------------------------------------------------------
    def finding_update_date(self, fl, update, update1, lats, size):
        """Retrieving the last update date"""
        s = ''
        update1[size:size+len(lats)] = s.join(update)

        return update1

#  ----------------------------------------------------------------------------
    def writing_variables(self):
        """Giving values to the variables"""
        #  Generation of the variables of the table

        lats, lons, juld, dac, wmo, update, prof, flag, size = self.reading_variables()
        self.latitudes[:] = lats[:size]
        self.longitudes[:] = lons[:size]
        self.juliandays[:] = juld[:size]
        self.dac[:] = dac[:size]
        self.wmo[:] = wmo[:size]
        self.update[:] = update[:size]
        self.prof[:] = prof[:size]
        self.flag[:] = flag[:size]

#  ----------------------------------------------------------------------------
    def creating_flags(self, f, flag, lats, size):
        """Creating the flags of the following profiles
         to know if they are OK or not"""
        #  flag = 0 : Seems to be OK
        #  flag = 101 : Variables missing
        #  flag = 102 : QC Position indicates wrong
        qc_pos = f.variables['POSITION_QC'][:]
        n_param = f.variables.has_key('PSAL') and f.variables.has_key('TEMP') and f.variables.has_key('PRES')
        for j, qc in enumerate(qc_pos):
            #  print(f.variables.has_key('PSAL'))
            if n_param:
                if qc == '1':
                    flag[size:size+j] = 0
                else:
                    flag[size:size+j] = 102
            else:
                flag[size:size+j] = 101
        return(flag)

#  ----------------------------------------------------------------------------
    def closing_file(self):
        """Closing the netCDF file"""
        #  Closing the netCDF file
        self.rootgrp.close()

#  ----------------------------------------------------------------------------

if __name__ == '__main__':

    #  Calling the constructor of the class
    __argodb__ = ArgoDB()
    __argodb__.closing_file()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
