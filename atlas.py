"""
Created on Mon Mar 12 13:10:24 2018
File creating the atlas of stats
from netCDF4  import Dataset
import numpy as np
import os
import argotools as argotools
import param as param
"""
from netCDF4  import Dataset
import numpy as np
import os
import argotools as argotools
import stats as stats
import param as param
import netCDF_form as ncform

dirstats = param.path_to_stats
diratlas = param.path_to_atlas

year = '2017'
month = '12'
day = '31'
date = [year, month, day]
# mode defines the values selected :
# R : Real time
# A : Adjusted Real Time
# D : Delayed time (Values verified)
mode = 'D'
typestat = 'zstd'
reso = 0.5
timeflag = 'annual'
nlon = 20
nlat = 15

# Will be used in the future to generate just the wanted variables
#  listvar = stats.var_choice

if typestat == 'zmean':
    listvar = ['NBbar', 'CTbar', 'SAbar', 'Ribar', 'BVF2bar']

elif typestat == 'zstd':
    listvar = ['NBstd', 'CTstd', 'SAstd', 'Ristd', 'BVF2std']

elif typestat == 'zdz':
    listvar = ['DZmean', 'DZstd', 'DZskew', 'EAPE']

else:
    raise ValueError('This typestat value does not exists')

dirstats = '%s/%g/%s/%s/%s' % (dirstats, reso, year, mode, typestat)

zref = argotools.zref

def ij2tile(i, j):
    return i + j*nlon

def gridindex2lonlat(ix, iy):
    
    lonmin_glo = -180.
    lonmax_glo = +180.
    latmin_glo = -80.
    latmax_glo = +80

    lon = lonmin_glo + ix*reso
    lat = latmin_glo + iy*reso
    return lon, lat

def get_glo_grid(reso):
    lonsize = []
    longlo = []
    latsize = []
    latglo = []
    for i in range(nlon):
        itile = i
        lat, lon = stats.grid_coordinate(itile, reso)
        lonsize.append(len(lon))
        longlo += list(lon)

    for j in range(nlat):
        itile = j*nlon
        lat, lon = stats.grid_coordinate(itile, reso)
        latsize.append(len(lat))
        latglo += list(lat)

    lonsize = np.asarray(lonsize)
    longlo = np.asarray(longlo)
    latsize = np.asarray(latsize)
    latglo = np.asarray(latglo)
    return lonsize, latsize, longlo, latglo
        
    
def atlas_filename(diratlas, reso, year, mode , typestat):

    if False: # many subfolders - short name for the atlas
        atlas_name = '%s_%g_annual' % (typestat, reso)
        ncfile = '%s/%s/%s/%s/%s/%s.nc' % (diratlas, reso, year, mode , typestat, atlas_name)

    else: # one folder - long name for the atlas
        ncfile ='%s/%s_%g_%s_%s.nc' % (diratlas, typestat, reso, year, mode)

    return ncfile


def glue_tiles(reso):
    """ Glue the stats tiles together into a global 3D atlas"""


    lonsize, latsize, lon, lat = get_glo_grid(reso)

    zref = argotools.zref
    nz = len(zref)

    ncfile = atlas_filename(diratlas, reso, year, mode , typestat)

    ncform.netCDF_dim_creation(ncfile, zref, len(lat), len(lon), mode, date)
    ncform.netCDF_var_creation(ncfile, stats.var_choice[typestat])

    lon2d, lat2d = np.meshgrid(lon, lat)

    nc = Dataset(ncfile, 'r+')
    nc.variables['zref'][:] = zref
    nc.variables['lon'][:, :] = lon2d
    nc.variables['lat'][:, :] = lat2d

    j0 = 0
    for j in range(nlat):
        i0 = 0
        for i in range(nlon):
            itile = ij2tile(i, j)
            i1 = i0+lonsize[i]
            j1 = j0+latsize[j]
            print('glue tile %3i: lon=%7.2f - lat=%7.2f' % (itile, lon[i0], lat[j0]))
            
            fname = stats.generate_filename(itile, typestat, reso, timeflag, date, mode)
            if os.path.isfile(fname):
                with Dataset(fname) as ncf:
                    for v in listvar:
                        z3d = ncf.variables[v][:, :, :]
                        for k in range(nz):
                            nc.variables[v][k, j0:j1, i0:i1] = z3d[k, :, :]
            else:
                for v in listvar:
                    for k in range(nz):
                        nc.variables[v][k, j0:j1, i0:i1] = 0.

            i0 = i1
        j0 = j1
    nc.close()
    print('Atlas %s has been generated successfully' % ncfile)



#  ----------------------------------------------------------------------------
if __name__ == '__main__':

    glue_tiles(reso)
