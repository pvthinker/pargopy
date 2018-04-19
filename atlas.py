from netCDF4  import Dataset
import numpy as np
import os
import argodb as argo
import argotools as argotools
import research_tools as res
import stats as stats

dirstats = '/net/libra/local/tmp/1/roullet/pargopy/stats'
diratlas = '/net/libra/local/tmp/1/roullet/pargopy/atlas'

dirstats = '/home2/datawork/therry/stats'
diratlas = '/home1/datawork/groullet/pargopy/atlas'

listvar = ['NBbar', 'CTbar', 'SAbar', 'Ribar']

reso = 0.25
nlon = 20
nlat = 15

dirstats = '%s/%g' % (dirstats, reso)
atlas_name = 'zmean_%g_annual' % reso

def ij2tile(i, j):
    return i + j*nlon



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
                        lonsize[i] = len(nc.variables['lon_deg'][0, :])
                        latsize[j]= len(nc.variables['lat_deg'][:, 0])

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
                    lon[i0:i1] = nc.variables['lon_deg'][0, :]
                    lat[j0:j1] = nc.variables['lat_deg'][:, 0]
                    zref = nc.variables['zref'][:]
            i0 = i1-1
        j0 = j1
    return lonsize, latsize, lon, lat, zref
    

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

    ncfile = '%s/%s.nc' % (diratlas, atlas_name)
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
lonsize, latsize, lon, lat, zref = get_glo_grid()
glue_tiles(np.arange(nlon),np.arange(nlat),lonsize, latsize, lon, lat, zref)
# glue_tiles(np.arange(5),np.arange(5),lonsize, latsize, lon, lat, zref)
