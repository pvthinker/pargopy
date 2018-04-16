"""Compute statistics on one tile"""
from netCDF4 import Dataset
import os
import numpy as np
import gsw as gsw
import time
from scipy import interpolate as ip
import general_tools as tools
import argodb as argo
import argotools as argotools
import tile as tiler
import param as param
import research_tools as research


path_to_stats = param.path_to_stats

zref = argotools.zref
daclist = argotools.daclist

def create_stat_file(itile, typestat, reso, timeflag):
    """Create statistics netcdf file"""

    rootgrp = Dataset('%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso, timeflag, itile), "w", format="NETCDF4")
    argodic = research.read_argo_filter(itile)
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso)
    nlat, nlon = np.shape(lon_deg)
    rootgrp.createDimension('depth', len(zref))
    rootgrp.createDimension('lat', nlat)
    rootgrp.createDimension('lon', nlon)

    if typestat == 'zmean':
        lon_deg = rootgrp.createVariable('lon_deg', 'f4', ('lat', 'lon'))
        lon_deg.long_name = 'Density'
        lon_deg.units = 'None'

        lat_deg = rootgrp.createVariable('lat_deg', 'f4', ('lat', 'lon'))
        lat_deg.long_name = 'Density'
        lat_deg.units = 'None'

        NBbar = rootgrp.createVariable('NBbar', 'f4', ('depth', 'lat', 'lon'))
        NBbar.long_name = 'Temperature'
        NBbar.units = 'Celsius'

        CTbar = rootgrp.createVariable('CTbar', 'f4', ('depth', 'lat', 'lon'))
        CTbar.long_name = 'Temperature'
        CTbar.units = 'Celsius'

        SAbar = rootgrp.createVariable('SAbar', 'f4', ('depth', 'lat', 'lon'))
        SAbar.long_name = 'Salinity'
        SAbar.units = '???'

        Ribar = rootgrp.createVariable('Ribar', 'f4', ('depth', 'lat', 'lon'))
        Ribar.long_name = 'Density'
        Ribar.units = 'None'

        zreference = rootgrp.createVariable('zref', 'f4', ('depth', ))
        zreference.long_name = 'Depth reference'
        zreference.units = 'Meter'

        rootgrp.close()

    elif typestat == 'zstd':

        lon_deg = rootgrp.createVariable('lon_deg', 'f4', ('lat', 'lon'))
        lon_deg.long_name = 'Density'
        lon_deg.units = 'None'

        lat_deg = rootgrp.createVariable('lat_deg', 'f4', ('lat', 'lon'))
        lat_deg.long_name = 'Density'
        lat_deg.units = 'None'

        NBstd = rootgrp.createVariable('NBstd', 'f4', ('depth', 'lat', 'lon'))
        NBstd.long_name = 'Temperature'
        NBstd.units = 'Celsius'

        CTstd = rootgrp.createVariable('CTstd', 'f4', ('depth', 'lat', 'lon'))
        CTstd.long_name = 'Density'
        CTstd.units = 'None'

        SAstd = rootgrp.createVariable('SAstd', 'f4', ('depth', 'lat', 'lon'))
        SAstd.long_name = 'Density'
        SAstd.units = 'None'

        Ristd = rootgrp.createVariable('Ristd', 'f4', ('depth', 'lat', 'lon'))
        Ristd.long_name = 'Density'
        Ristd.units = 'None'

        v = rootgrp.createVariable('DZstd', 'f4', ('depth', 'lat', 'lon'))
        v.long_name = 'Isopycnal displacement'
        v.units = 'm'

        v = rootgrp.createVariable('EAPE', 'f4', ('depth', 'lat', 'lon'))
        v.long_name = 'Eddy available potential energy'
        v.units = 'J.m^-3'

        rootgrp.close()

    # nom de fichier:
    # typestat_reso_timeflag_tileidx.nc
    #
    # where
    # typestat = 'zmean', 'zstd' etc
    # reso = 0.5, 1.0, 0.25 etc
    # timeflag = 'annual', 'DJF', 'MAM', 'JJA', 'SON'
    # tileidx = 0 ... 299

    # en fonction du type de stat on ne sauve pas les memes variables...

    # pour ecrire un attribut global_clim
    # nc.setncattr('lonmin', -32.)if typestat == 'zmean':


def write_stat_file(itile, typestat, reso_deg, timeflag):
    """Write statistics into a netcdf file"""
    # idem, depend du type de stat
    filename = '%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso_deg, timeflag, itile)
    if (os.path.isfile(filename)):
        print('filename existe')
        f = Dataset(filename, "r+", format="NETCDF4")
        # idem, depend du type de stat
        if typestat == 'zmean':
            lon_deg, lat_deg, NBbar, CTbar, SAbar, Ribar = compute_mean_at_zref(itile, reso_deg)
            f.variables['CTbar'][:, :, :] = CTbar
            f.variables['SAbar'][:, :, :] = SAbar
            f.variables['Ribar'][:, :, :] = Ribar
            f.variables['NBbar'][:, :, :] = NBbar
            f.variables['lat_deg'][:, :] = lat_deg
            f.variables['lon_deg'][:, :] = lon_deg
            f.close()
        elif typestat == 'zstd':
            lon_deg, lat_deg, CTstd, SAstd, DZstd, Ristd, EAPE, NBstd = compute_std_at_zref(itile, reso_deg, timeflag)
            f.variables['NBstd'][:, :, :] = NBstd
            f.variables['CTstd'][:, :, :] = CTstd
            f.variables['SAstd'][:, :, :] = SAstd
            f.variables['DZstd'][:, :, :] = DZstd
            f.variables['Ristd'][:, :, :] = Ristd
            f.variables['EAPE'][:, :, :] = EAPE
            f.variables['lat_deg'][:, :] = lat_deg
            f.variables['lon_deg'][:, :] = lon_deg
            f.close()


def read_stat_file(typestat, itile, reso, timeflag):
    """Read statistics into a netcdf file"""
    filename = '%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso, timeflag, itile)
    if (os.path.isfile(filename)):
        f = Dataset(filename, "r", format="NETCDF4")
        # idem, depend du type de stat
        if typestat == 'zmean':
            CTbar = f.variables['CTbar'][:, :, :]
            SAbar = f.variables['SAbar'][:, :, :]
            Ribar = f.variables['Ribar'][:, :, :]
            NBbar = f.variables['NBbar'][:, :, :]
            lat_deg = f.variables['lat_deg'][:, :]
            lon_deg = f.variables['lon_deg'][:, :]
            f.close()
            return lon_deg, lat_deg, NBbar, CTbar, SAbar, Ribar
        elif typestat == 'zstd':
            NBstd = f.variables['CTstd'][:, :, :]
            CTstd = f.variables['CTstd'][:, :, :]
            SAstd = f.variables['SAstd'][:, :, :]
            Ristd = f.variables['Ristd'][:, :, :]
            lat_deg = f.variables['lat_deg'][:, :]
            lon_deg = f.variables['lon_deg'][:, :]
            f.close()
            return lon_deg, lat_deg, CTstd, SAstd, Ristd, NBstd


def define_grid(minlon, maxlon, minlat, maxlat, reso_deg):
    """ setup the grid coordinates (in degrees)
    coordinates are round multiples of reso_deg
    reso_deg sets the grid resolution, typically 0.5deg""" 

    minlon = np.ceil(minlon/reso_deg)*reso_deg
    maxlon = np.floor(maxlon/reso_deg)*reso_deg

    minlat = np.ceil(minlat/reso_deg)*reso_deg
    maxlat = np.floor(maxlat/reso_deg)*reso_deg

    lon1D = np.arange(minlon, maxlon+reso_deg, reso_deg)
    lat1D = np.arange(minlat, maxlat+reso_deg, reso_deg)

    lon, lat = np.meshgrid(lon1D, lat1D)

    return lon, lat


def compute_mean_at_zref(itile, reso_deg):
    """Compute the mean at depths zref"""
    # argodb = research.read_argo_filter(itile)  # argo.read_argodb()
    tile = tiler.read_tile(itile)
    # output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    # nbprof = output['IPROF']
    CT, SA, RI, lat, lon = tile['CT'], tile['SA'], tile['RHO'], tile['LATITUDE'], tile['LONGITUDE']
    argodic = research.read_argo_filter(itile)
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso_deg)

    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(zref)
    nbprof = len(CT)

    # gridded arrays of CT, SA et RI means
    NBbar = np.zeros((nz, nlat, nlon))
    CTbar = np.zeros((nz, nlat, nlon))
    SAbar = np.zeros((nz, nlat, nlon))
    RIbar = np.zeros((nz, nlat, nlon))
    for k in range(len(nbprof)):
        print('%4i/%i' % (k, len(nbprof)))
        # todo: weigh in time using juld,
        # e.g. only winter statistics
        time_weight = 1.

        xlon_rad = np.deg2rad(lon[k])
        xlat_rad = np.deg2rad(lat[k])
        weight = tools.compute_weight(lon_rad, lat_rad,
                                      xlon_rad, xlat_rad,
                                      reso_rad)
        weight *= time_weight

        for l in range(nz):
            if np.isnan(CT[k, l]) or np.isnan(SA[k, l]):
                pass
            else:
                NBbar[l, :, :] += weight
                CTbar[l, :, :] += weight*CT[k, l]
                SAbar[l, :, :] += weight*SA[k, l]
                RIbar[l, :, :] += weight*RI[k, l]

    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    coef = 1./NBbar
    coef[NBbar < 1] = np.NaN

    CTbar *= coef
    SAbar *= coef
    RIbar *= coef

    return lon_deg, lat_deg, NBbar, CTbar, SAbar, RIbar


def compute_std_at_zref(itile, reso_deg, timeflag, verbose=False):
    """Compute the standard deviations at depths zref"""

    # gridded arrays of CT, SA variances
    lon_deg, lat_deg, NBbar, CAbar, SAbar, RHObar = read_stat_file('zmean', itile, reso_deg, timeflag) # read it from the file
    argodb = research.read_argo_filter(itile)  # argo.read_argodb()
    tile = tiler.read_tile(itile)
    # output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    # iprof = output['IPROF']
    CT, SA, RI, lat, lon = tile['CT'], tile['SA'], tile['RHO'], tile['LATITUDE'], tile['LONGITUDE']

    nbprof = len(CT)
    
    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(lon_deg)

    # RI is rho in situ

    nz = len(zref)

    xlon_rad = np.deg2rad(lon)
    xlat_rad = np.deg2rad(lat)

    NBstd = np.zeros((nz, nlat, nlon))
    CTstd = np.zeros((nz, nlat, nlon))
    SAstd = np.zeros((nz, nlat, nlon))
    DZstd = np.zeros((nz, nlat, nlon))
    Ristd = np.zeros((nz, nlat, nlon))
    EAPE = np.zeros((nz, nlat, nlon))

    wmin = 5e-3 # minimum weight below which a profile is drop

    # double loop on each grid point (instead of a loop on each profile)
    for j in range(nlat):
        for i in range(nlon):
            if verbose:
                print('%i/%i - %i/%i' % (i, nlon, j, nlat))
            time_weight = 1.

            # do all profiles
            weight = tools.compute_weight(lon_rad[j, i], lat_rad[j, i],
                                            xlon_rad, xlat_rad,
                                            reso_rad)
            weight *= time_weight
            # print(np.shape(RHObar))
            interpolator = ip.interp1d(RHObar[:,j,i], zref, bounds_error=False)

            p = gsw.p_from_z(-zref, lat[j])
            g = gsw.grav(lat[j], p)
            cs = gsw.sound_speed(SAbar[:, j, i], CAbar[:, j, i], p)
            # the copy is to make data contiguous in memory
            rho0 = RHObar[:, j, i].copy()

            for l in range(nbprof):
                if weight[l] < wmin:
                    pass
                else:
                    # the vertical isopycnal displacement dz
                    # is computed by first
                    # interpolating RI(z) on the reference profile
                    # z0(rho) [inverse mapping of rho0(z)]
                    # we get dzstar
                    # then by adding a correction due to
                    # compressibility
                    dzstar = interpolator(RI[l, :]) - zref
                    drho = RI[l, :] - rho0
                    # add correction
                    dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
                    eape = 0.5*dz*drho

                    for k in range(nz):
                        if np.isnan(dz[k]) or np.isnan(drho[k]):
                            pass
                        else:
                            NBstd[k, j, i] += weight[l]
                            CTstd[k, j, i] += weight[l]*CT[l, k]**2
                            SAstd[k, j, i] += weight[l]*SA[l, k]**2
                            DZstd[k, j, i] += weight[l]*dz[k]**2
                            Ristd[k, j, i] += weight[l]*drho[k]**2
                            EAPE[k, j, i] += weight[l]*dz[k]*drho[k]
    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    # std = sqrt( 1/(n-1) sum_i (x_i -xbar)^2)
    #     = sqrt( 1/(n-1) sum_i x_i^2 - n/(n-1)*xbar^2)
    coef = 1./(NBstd-1)
    coef[NBstd < 2] = 0.

    CTstd = np.sqrt(coef*(CTstd-NBstd*CAbar))
    SAstd = np.sqrt(coef*(SAstd-NBstd*SAbar))
    DZstd = np.sqrt(coef*DZstd)
    Ristd = np.sqrt(coef*Ristd)
    EAPE = 0.5*coef*EAPE

    return lon_deg, lat_deg, CTstd, SAstd, DZstd, Ristd, EAPE, NBstd


def main(itile, typestat, reso, timeflag):
    """Main function of stats.py"""    
    create_stat_file(itile, typestat, reso, timeflag)
    write_stat_file(itile, typestat, reso, timeflag)
    #  create_stat_file(68, 'zmean', 0.5, 'annual')
    #  write_stat_file(68, 'zmean', 0.5, 'annual')


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    main(50, 'zstd', 0.5, 'annual')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
