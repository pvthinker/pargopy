# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:14:57 2018

@author: herry
"""

from netCDF4 import Dataset
import numpy as np
import os
import math
import time
import research_tools as research
import interpolation_tools as interpolation
import param as param
tmps1 = time.time()


class Tile:
    """Class used to classify the rootgrps with their latitude/longitude
        into boxes"""

    def __init__(self, tile_idx):
        """Initialization of the tiles"""
        # List of the datacenters
        self.daclist=['aoml','bodc','coriolis','csio',
                  'csiro','incois','jma','kma',
                  'kordi','meds','nmdis']
        #  Path to ARGO
        self.path, self.path_to_tiles, self.path_to_argotiles = param.data_tile()
        # Normalized depth
        self.zref = np.array([0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
                      110., 120., 130., 140., 150., 160., 170., 180., 190.,
                      200., 220., 240., 260., 280, 300., 320., 340., 360.,
                      380., 400., 450., 500., 550., 600., 650., 700., 750.,
                      800., 850., 900., 950., 1000., 1050., 1100., 1150.,
                      1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550.,
                      1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950.,
                      2000.])
        #  Max number of temperature, salinity and pressure values
        self.len_values = len(self.zref)
        # Index of the tile
        self.tile_idx = tile_idx

#  ----------------------------------------------------------------------------
    def tile_creation(self):
        #  Creating the following netCDF chunk files
        self.rootgrp = Dataset(self.path_to_tiles % self.tile_idx, "w", format="NETCDF4")
        self.variables_creation()

#  ----------------------------------------------------------------------------
    def variables_creation(self):
        """Creating the dimensions of the variables"""
        #  Generation of the dimension of the table
        self.rootgrp.createDimension('nb_prof', None)
        self.rootgrp.createDimension('lev', len(self.zref))

        #  Creation of the variables of the table (T°, Salinity, Pressure)

        self.CT = self.rootgrp.createVariable('CT', 'f4', ('nb_prof', 'lev'))
        self.CT.long_name = 'Temperature'
        self.CT.units = 'Celsius'

        self.SA = self.rootgrp.createVariable('SA', 'f4', ('nb_prof', 'lev'))
        self.SA.long_name = 'Salinity'
        self.SA.units = '???'

        self.Rho = self.rootgrp.createVariable('Rho', 'f4', ('nb_prof', 'lev'))
        self.Rho.long_name = 'Density'
        self.Rho.units = 'None'

        self.juld = self.rootgrp.createVariable("juld", "f4", ('nb_prof'))
        self.juld.long_name = 'julian day'
        self.juld.units = 'degrees'

        self.latitudes = self.rootgrp.createVariable("lat", "f4", ('nb_prof'))
        self.latitudes.long_name = 'latitude'
        self.latitudes.units = 'degrees'

        self.longitudes = self.rootgrp.createVariable("lon", "f4", ('nb_prof'))
        self.longitudes.long_name = 'longitude'
        self.longitudes.units = 'degrees'

        self.prof = self.rootgrp.createVariable('prof', 'i2', ('nb_prof', ))
        self.prof.long_name = 'index of the profile in the profiles dir'
        self.prof.units = 'dimensionless'

        self.dac = self.rootgrp.createVariable('dac', 'b', ('nb_prof', ))
        self.dac.long_name = 'index of the dac'
        self.dac.units = 'dimensionless'

        self.wmo = self.rootgrp.createVariable('wmo', 'i', ('nb_prof', ))
        self.wmo.long_name = 'wmoid of the wmo'
        self.wmo.units = 'dimensionless'

        self.zreference = self.rootgrp.createVariable('zref', 'f4', ('lev', ))
        self.zreference.long_name = 'Depth reference'
        self.zreference.units = 'Meter'

        return (self.CT, self.SA, self.Rho, self.juld, self.latitudes, self.longitudes, self.dac, self.wmo, self.prof, self.zreference)

#  ----------------------------------------------------------------------------
    def reading_variables(self, i):
        """Salvation of the profiles datas"""
        #  Reading the netCDF file fom ARGO
        if (os.path.isfile(self.path_to_argotiles % i)):
            f = Dataset(self.path_to_argotiles % i, "r", format="NETCDF4")
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            juld = f.variables['juld'][:]
            prof = f.variables['prof'][:]
            dac = f.variables['dac'][:]
            wmo = f.variables['wmo'][:]
            latmin = f.latmin
            latmax = f.latmax
            lonmin = f.lonmin
            lonmax = f.lonmax
            line_number = f.dimensions['line_number'].size

            f.close()

        return (prof, dac, wmo, lats, lons, juld, line_number, latmin, latmax, lonmin, lonmax)

#  ----------------------------------------------------------------------------
    def reading_ARGO(self, dac, wmo):
        """Reading the netCDF file from ARGO corresponding to the dac/wmo founded"""

        filename = '%s/%s/%s/%s_prof.nc' % (self.path, self.daclist[dac], wmo, wmo)
        print(filename)
        if (os.path.isfile(filename)):
            f = Dataset(filename, "r", format="NETCDF4")
            temp = f.variables['TEMP'][:, :]
            temp_qc = f.variables['TEMP_QC'][:, :]
            sal = f.variables['PSAL'][:, :]
            sal_qc = f.variables['PSAL_QC'][:, :]
            pres = f.variables['PRES'][:, :]
            pres_qc = f.variables['PRES_QC'][:, :]
            lat = f.variables['LATITUDE'][:]
            lon = f.variables['LONGITUDE'][:]
            f.close()

        return temp, sal, pres, temp_qc, sal_qc, pres_qc, lat, lon

#  ----------------------------------------------------------------------------
    def reading_tile(self):
        """Reading the netCDF file describing the tiles"""
        filename = self.path_to_argotiles % self.tile_idx
        if (os.path.isfile(filename)):
            f = Dataset(filename, "r", format="NETCDF4")
            CT = f.variables['CT'][:, :]
            SA = f.variables['SA'][:, :]
            Rho = f.variables['Rho'][:, :]
            prof = f.variables['prof'][:]
            lat = f.variables['lat'][:]
            lon = f.variables['lon'][:]
            dac = f.variables['dac'][:]
            wmo = f.variables['wmo'][:]
            f.close()

        return CT, SA, prof, lat, lon, dac, wmo

#  ----------------------------------------------------------------------------
    def writing_variables(self, Ti, Si, Ri, lats1, lons1, i):
        """Giving values to the variables"""
        #  Generation of the values into the tiles
        prof, dac, wmo, lats, lons, juld, line_number, latmin, latmax, lonmin, lonmax = self.reading_variables(i)
        print(lats1)
        print(lats)
        self.CT[:, :] = Ti
        self.SA[:, :] = Si
        self.Rho[:, :] = Ri
        self.latitudes[:] = lats1
        self.longitudes[:] = lons1
        self.dac[:] = dac
        self.wmo[:] = wmo
        self.prof[:] = prof
        self.zreference[:] = self.zref
        self.juld[:] = juld
        self.rootgrp.setncattr('latmin', latmin)
        self.rootgrp.setncattr('latmax', latmax)
        self.rootgrp.setncattr('lonmin', lonmin)
        self.rootgrp.setncattr('lonmax', lonmax)

#  ----------------------------------------------------------------------------
    def interpolate_all_profiles(self):
        """Giving values to the variables"""
        #  Generation of the values into the tiles

        #  Reading the netCDF file fom ARGO
        #  prof correponds to the index of the profile
        prof, dac, wmo, lats, lons, juld, line_number, latmin, latmax, lonmin, lonmax = self.reading_variables(self.tile_idx)
        len_argotile = len(prof)  # line_number
        lats1 = np.zeros((len_argotile,)) + np.NaN
        lons1 = np.zeros((len_argotile,)) + np.NaN
        juld1 = np.zeros((len_argotile,)) + np.NaN
        prof1 = np.zeros((len_argotile,)) + np.NaN
        temp1 = np.zeros((len_argotile, self.len_values)) + np.NaN
        sal1 = np.zeros((len_argotile, self.len_values)) + np.NaN
        rho1 = np.zeros((len_argotile, self.len_values)) + np.NaN
        wmoid = list(set(wmo))
        k = 0
        for j, xwmoid in enumerate(wmoid):
            dacid = research.retrieve_dac_from_wmoid(dac, wmo, wmoid[j])
            prof1 = research.retrieve_prof_from_wmoid(prof, wmo, wmoid[j])
            temp, sal, pres, temp_qc, sal_qc, pres_qc, lat, lon = self.reading_ARGO(dacid, wmoid[j])
            for i, idx in enumerate(prof1):
                Ti, Si, Ri, CT, SA, z, ierr = interpolation.raw_to_interpolate(temp[idx], sal[idx], pres[idx], temp_qc[idx], sal_qc[idx], pres_qc[idx], lon[idx], lat[idx], self.zref)
                if ierr != 0:
                    print('ierr != 0')
                else:
                    #  lats1[k] = lat[i]
                    #  lons1[k] = lon
                    #  juld1[k] = juld
                    temp1[k, :len(Ti)] = Ti
                    sal1[k, :len(Si)] = Si
                    rho1[k, :len(Ri)] = Ri
                    k += 1
        #  print(len(temp1))
        self.writing_variables(temp1, sal1, rho1, lats1, lons1, self.tile_idx)

#  ----------------------------------------------------------------------------
    def closing_file(self):
        """Closing the netCDF file"""
        #  Closing the netCDF file
        self.rootgrp.close()
#  ----------------------------------------------------------------------------
if __name__ == '__main__':
#==============================================================================
#     __tile__ = Tile(51)
#     __tile__.tile_creation()
#     __tile__.interpolate_all_profiles()
#     __tile__.closing_file()
#==============================================================================
#==============================================================================
#     __tile__ = Tile(226)
#     __tile__.tile_creation()
#     __tile__.interpolate_all_profiles()
#     __tile__.closing_file()
#==============================================================================
#==============================================================================
#     __tile__ = Tile(140)
#     __tile__.tile_creation()
#     __tile__.interpolate_all_profiles()
#     __tile__.closing_file()
#==============================================================================
#==============================================================================
#     __tile__ = Tile(213)
#     __tile__.tile_creation()
#     __tile__.interpolate_all_profiles()
#     __tile__.closing_file()
#==============================================================================
    __tile__ = Tile(172)
    __tile__.tile_creation()
    __tile__.interpolate_all_profiles()
    __tile__.closing_file()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
