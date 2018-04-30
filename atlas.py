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
import param as param

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

if typestat == 'zmean':
    listvar = ['NBbar', 'CTbar', 'SAbar', 'Ribar', 'BVF2bar']
elif typestat == 'zstd':
    listvar = ['NBstd', 'CTstd', 'SAstd', 'Ristd', 'BVF2std', 'DZmean', 'DZstd', 'DZskew', 'EAPE']
else:
    raise ValueError('This typestat value does not exists')

dirstats = '%s/%g/%s/%s/%s' % (dirstats, reso, year, mode, typestat)
atlas_name = '%s_%g_annual' % (typestat, reso)

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

def get_glo_grid():
    lonsize = np.zeros((nlon,), dtype=int)
    latsize = np.zeros((nlat,), dtype=int)

    lonmin_glo = -180.
    lonmax_glo = +180.
    latmin_glo = -80.
    latmax_glo = +80

    nlon_glo = int((lonmax_glo-lonmin_glo)/reso)+1
    nlat_glo = int((latmax_glo-latmin_glo)/reso)+1
    lon = np.zeros((nlon_glo,))
    lat = np.zeros((nlat_glo,))

    print('global grid has %i x %i points' % (nlon_glo, nlat_glo))
    
    for j in range(nlat):
        for i in range(nlon):
            if (latsize[j]==0) or (lonsize[i]==0):
                itile = ij2tile(i, j)
                ncfile = '%s/%s_%03i.nc' % (dirstats, atlas_name, itile)
                if os.path.isfile(ncfile):
                    with Dataset(ncfile, 'r') as nc:
                        lonsize[i] = len(nc.variables['lon'][0, :])
                        latsize[j]= len(nc.variables['lat'][:, 0])

    j0 = 0
    for j in range(nlat):
        j1 = j0+latsize[j]
        i0 = 0
        for i in range(nlon):
            i1 = i0+lonsize[i]
            
            itile = ij2tile(i, j)
            ncfile = '%s/%s_%03i.nc' % (dirstats, atlas_name, itile)
            if os.path.isfile(ncfile):
                print(i,j,itile)
                with Dataset(ncfile, 'r') as nc:
                    lon[i0:i1] = nc.variables['lon'][0, :]
                    lat[j0:j1] = nc.variables['lat'][:, 0]
                    zref = nc.variables['zref'][:]
            i0 = i1-1
        j0 = j1
    return lonsize, latsize, lon, lat
    

def glue_tiles(iloc, jloc,lonsize, latsize, lon, lat, zref):
    # lon = []
    # j = jloc[0]
    # for i in iloc:
    #     itile = ij2tile(i, j)
    #     ncfile = '%s/%s_%03i.nc' % (dirstats, atlas_name, itile)
    #     with Dataset(ncfile, 'r') as nc:
    #         lon += list(nc.variables['lon_deg'][0, :])

    # lat = []
    # i = iloc[0]
    # for j in jloc:
    #     itile = ij2tile(i, j)
    #     ncfile = '%s/%s_%03i.nc' % (dirstats, atlas_name, itile)
    #     with Dataset(ncfile, 'r') as nc:
    #         lat += list(nc.variables['lat_deg'][:, 0])
            
    # with Dataset(ncfile, 'r') as nc:
    #     zref = nc.variables['zref'][:]

    print(lon)
    print(lat)
    print(zref)
    
    nz = len(zref)

    ncfile = '%s/%s/%s/%s/%s/%s.nc' % (diratlas, reso, year, mode , typestat, atlas_name)
    nc = Dataset(ncfile, 'w', format='NETCDF4')
    nc.createDimension('zref', nz)
    nc.createDimension('lat', len(lat))       
    nc.createDimension('lon', len(lon))
    var = nc.createVariable('zref', 'f', ('zref'))
    var.long_name = 'zref'
    var = nc.createVariable('lon', 'f', ('lon'))
    var.long_name = 'lon'
    var = nc.createVariable('lat', 'f', ('lat'))
    var.long_name = 'lat'
    for v in listvar:
        var = nc.createVariable(v,'f',('zref','lat','lon'))
        var.long_name = v
    nc.close()

    nc = Dataset(ncfile, 'r+')
    nc.variables['zref'][:] = zref
    nc.variables['lon'][:] = lon
    nc.variables['lat'][:] = lat

    j0 = 0
    for j in jloc:
        i0 = 0
        for i in iloc:
            itile = ij2tile(i, j)
            i1 = i0+lonsize[i]  # len(ncf.dimensions['lon'])
            j1 = j0+latsize[j]  # len(ncf.dimensions['lat'])
            print('tile %i : lon=%g - lat=%g' % (itile,lon[i0], lat[j0]))
            fname = '%s/%s_%03i.nc' % (dirstats, atlas_name, itile)
            if os.path.isfile(fname):
                with Dataset(fname) as ncf:
                    for v in listvar:
                        z3d = ncf.variables[v][:, : :]
                        for k in range(nz):
                            nc.variables[v][k, j0:j1, i0:i1] = z3d[k, :, :]
            else:
                pass
                # for v in listvar:
                #     for k in range(nz):
                #         nc.variables[v][k, j0:j1, i0:i1] = 0.

            i0 = i1-1
        j0 = j1
    nc.close()

# glue_tiles(np.arange(10,13),np.arange(2,4))


# glue_tiles(np.arange(5),np.arange(5),lonsize, latsize, lon, lat, zref)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    lonsize, latsize, lon, lat = get_glo_grid()
    glue_tiles(np.arange(nlon),np.arange(nlat),lonsize, latsize, lon, lat, zref)