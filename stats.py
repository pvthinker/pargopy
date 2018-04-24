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
    filename = '%s/%s/%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso, typestat, reso, timeflag, itile)
    rootgrp = Dataset(filename, "w", format="NETCDF4")
    argodic = argotools.read_argo_filter(itile)
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso)
    nlat, nlon = np.shape(lon_deg)
    rootgrp.createDimension('zref', len(zref))
    rootgrp.createDimension('lat', nlat)
    rootgrp.createDimension('lon', nlon)

    if typestat == 'zmean':
        lon_deg = rootgrp.createVariable('lon', 'f4', ('lat', 'lon'))
        lon_deg.long_name = 'Longitude in degrees'
        lon_deg.units = 'Degrees'

        lat_deg = rootgrp.createVariable('lat', 'f4', ('lat', 'lon'))
        lat_deg.long_name = 'Latitude in degrees'
        lat_deg.units = 'Degrees'

        NBbar = rootgrp.createVariable('NBbar', 'f4', ('zref', 'lat', 'lon'))
        NBbar.long_name = 'Weight'
        NBbar.units = 'None'

        CTbar = rootgrp.createVariable('CTbar', 'f4', ('zref', 'lat', 'lon'))
        CTbar.long_name = 'Temperature'
        CTbar.units = 'Celsius'

        SAbar = rootgrp.createVariable('SAbar', 'f4', ('zref', 'lat', 'lon'))
        SAbar.long_name = 'Salinity'
        SAbar.units = '???'

        Ribar = rootgrp.createVariable('Ribar', 'f4', ('zref', 'lat', 'lon'))
        Ribar.long_name = 'Density'
        Ribar.units = '???'

        BVF2bar = rootgrp.createVariable('BVF2bar', 'f4', ('zref', 'lat', 'lon'))
        BVF2bar.long_name = 'Brunt Vaisala frequency squared'
        BVF2bar.units = 's^-2'

        zreference = rootgrp.createVariable('zref', 'f4', ('zref', ))
        zreference.long_name = 'Depth reference'
        zreference.units = 'Meter'

        rootgrp.close()

    elif typestat == 'zstd':

        zreference = rootgrp.createVariable('zref', 'f4', ('zref', ))
        zreference.long_name = 'Depth reference'
        zreference.units = 'Meter'

        lon_deg = rootgrp.createVariable('lon', 'f4', ('lat', 'lon'))
        lon_deg.long_name = 'Longitude in degrees'
        lon_deg.units = 'Degrees'

        lat_deg = rootgrp.createVariable('lat', 'f4', ('lat', 'lon'))
        lat_deg.long_name = 'Latitude in degrees'
        lat_deg.units = 'Degrees'

        NBstd = rootgrp.createVariable('NBstd', 'f4', ('zref', 'lat', 'lon'))
        NBstd.long_name = 'Weight'
        NBstd.units = 'None'

        CTstd = rootgrp.createVariable('CTstd', 'f4', ('zref', 'lat', 'lon'))
        CTstd.long_name = 'Temperature'
        CTstd.units = 'Degree Celsius'

        SAstd = rootgrp.createVariable('SAstd', 'f4', ('zref', 'lat', 'lon'))
        SAstd.long_name = 'Salinity'
        SAstd.units = 'g.kg^-1'

        Ristd = rootgrp.createVariable('Ristd', 'f4', ('zref', 'lat', 'lon'))
        Ristd.long_name = 'Density'
        Ristd.units = 'kg.m^-3'

        BVF2std = rootgrp.createVariable('BVF2std', 'f4', ('zref', 'lat', 'lon'))
        BVF2std.long_name = 'std of Brunt Vaisala frequency squared'
        BVF2std.units = 's^-2'

        v = rootgrp.createVariable('DZmean', 'f4', ('depth', 'lat', 'lon'))
        v.long_name = 'Mean isopycnal displacement'
        v.units = 'm'

        v = rootgrp.createVariable('DZstd', 'f4', ('depth', 'lat', 'lon'))
        v.long_name = 'Std isopycnal displacement'
        v.units = 'm'

        v = rootgrp.createVariable('DZskew', 'f4', ('depth', 'lat', 'lon'))
        v.long_name = 'Skewness isopycnal displacement'
        v.units = 'none'

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
    filename = '%s/%s/%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso_deg, typestat, reso_deg, timeflag, itile)
    if (os.path.isfile(filename)):
        print('filename existe')
        f = Dataset(filename, "r+", format="NETCDF4")
        # idem, depend du type de stat
        if typestat == 'zmean':
            lon_deg, lat_deg, NBbar, CTbar, SAbar, Ribar, BVF2bar = compute_mean_at_zref(itile, reso_deg)
            f.variables['CTbar'][:, :, :] = CTbar
            f.variables['SAbar'][:, :, :] = SAbar
            f.variables['Ribar'][:, :, :] = Ribar
            f.variables['BVF2bar'][:, :, :] = BVF2bar
            f.variables['NBbar'][:, :, :] = NBbar
            f.variables['lat'][:, :] = lat_deg
            f.variables['lon'][:, :] = lon_deg
            f.variables['zref'][:] = zref
            f.close()
        elif typestat == 'zstd':
            lon_deg, lat_deg, CTstd, SAstd, BVF2std, DZmean, DZstd, DZskew, Ristd, EAPE, NBstd = compute_std_at_zref(itile, reso_deg, timeflag)
            f.variables['NBstd'][:, :, :] = NBstd
            f.variables['CTstd'][:, :, :] = CTstd
            f.variables['SAstd'][:, :, :] = SAstd
            f.variables['BVF2std'][:, :, :] = BVF2std
            f.variables['DZmean'][:, :, :] = DZmean
            f.variables['DZstd'][:, :, :] = DZstd
            f.variables['DZskew'][:, :, :] = DZskew
            f.variables['Ristd'][:, :, :] = Ristd
            f.variables['EAPE'][:, :, :] = EAPE
            f.variables['lat'][:, :] = lat_deg
            f.variables['lon'][:, :] = lon_deg
            f.variables['zref'][:] = zref
            f.close()


def read_stat_file(typestat, itile, reso, timeflag):
    """Read statistics into a netcdf file"""
    filename = '{0}/{1}/{2}/{1}_{2}_{3}_{4:03}.nc'.format(path_to_stats, typestat, reso, timeflag, itile)
    print('read stat file : %s' % filename)
    # filename = '%s/%s/%s/%s_%s_%s_%003i.nc' % (path_to_stats, typestat, reso, typestat, reso, timeflag, itile)

    if (os.path.isfile(filename)):
        f = Dataset(filename, "r", format="NETCDF4")
        # idem, depend du type de stat
        if typestat == 'zmean':
            CTbar = f.variables['CTbar'][:, :, :]
            SAbar = f.variables['SAbar'][:, :, :]
            Ribar = f.variables['Ribar'][:, :, :]
            #BVF2bar = f.variables['BVF2bar'][:, :, :]
            NBbar = f.variables['NBbar'][:, :, :]
            lat_deg = f.variables['lat'][:, :]
            lon_deg = f.variables['lon'][:, :]
            f.close()
            return lon_deg, lat_deg, NBbar, CTbar, SAbar, Ribar#, BVF2bar
        elif typestat == 'zstd':
            NBstd = f.variables['CTstd'][:, :, :]
            CTstd = f.variables['CTstd'][:, :, :]
            SAstd = f.variables['SAstd'][:, :, :]
            Ristd = f.variables['Ristd'][:, :, :]
            # todo: BVF2std
            lat_deg = f.variables['lat'][:, :]
            lon_deg = f.variables['lon'][:, :]
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
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']

    argodic = argotools.read_argo_filter(itile)
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
    BVF2bar = np.zeros((nz, nlat, nlon))
    for k in range(nbprof):
        print('%4i/%i' % (k, nbprof))
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
                BVF2bar[l, :, :] += weight*BVF2[k, l]

    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    coef = 1./NBbar
    coef[NBbar < 1] = np.NaN

    CTbar *= coef
    SAbar *= coef
    RIbar *= coef
    BVF2bar *= coef

    return lon_deg, lat_deg, NBbar, CTbar, SAbar, RIbar, BVF2bar


def compute_std_at_zref(itile, reso_deg, timeflag, verbose=False):
    """Compute the standard deviations at depths zref"""

    # gridded arrays of CT, SA variances
    lon_deg, lat_deg, NBbar, CAbar, SAbar, RHObar = read_stat_file('zmean', itile, reso_deg, timeflag) # read it from the file
    argodb = argotools.read_argodb()
    tile = tiler.read_tile(itile)
    # output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    # iprof = output['IPROF']
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']

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
    BVF2std = np.zeros((nz, nlat, nlon))
    DZmean = np.zeros((nz, nlat, nlon))
    DZstd = np.zeros((nz, nlat, nlon))
    DZskew = np.zeros((nz, nlat, nlon))
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
            rho0 = RHObar[:, j, i].copy()

            drho = RI - RHObar[:, j, i]
            dCT = CT-CAbar[:, j, i]
            dSA = SA-SAbar[:, j, i]
            zrho = interpolator(RI)
            dzstar = zrho-zref
            dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
            eape = 0.5*dz*drho

            weight = weight[:, np.newaxis] + np.zeros_like(zref)
            weight[np.where(np.isnan(dz) | np.isnan(drho))] = 0.

            def average(field):
                return np.nansum(weight*field, axis=0)
            
            NBstd[:, j, i] = average(1.)
            CTstd[:, j, i] = average(dCT**2)
            SAstd[:, j, i] = average(dSA**2)
            DZmean[:, j, i] = average(dz)
            DZstd[:, j, i] = average(dz**2)
            DZskew[:, j, i] = average(dz**3)
            Ristd[:, j, i] = average(drho**2)
            EAPE[:, j, i] = average(dz*drho)
            # NBstd[:, j, i] = np.nansum(weight, axis=0)
            # CTstd[:, j, i] = np.nansum(weight*dCT**2, axis=0)
            # SAstd[:, j, i] = np.nansum(weight*dSA**2, axis=0)
            # DZmean[:, j, i] = np.nansum(weight*dz, axis=0)
            # DZstd[:, j, i] = np.nansum(weight*dz**2, axis=0)
            # DZskew[:, j, i] = np.nansum(weight*dz**3, axis=0)
            # Ristd[:, j, i] = np.nansum(weight*drho**2, axis=0)
            # EAPE[:, j, i] = np.nansum(weight*dz*drho, axis=0)

    # normalize with the number of profiles (fractional
    # because NBbar is fractionnal)
    # std = sqrt( 1/(n-1) sum_i (x_i -xbar)^2)
    #     = sqrt( 1/(n-1) sum_i x_i^2 - n/(n-1)*xbar^2)
    coef = 1./(NBstd-1)
    coef[NBstd < 2] = np.nan

    CTstd = np.sqrt(coef*CTstd)
    SAstd = np.sqrt(coef*SAstd)
    DZmean *= coef
    DZstd = np.sqrt(coef*DZstd)
    # skew = E( ((X-mu)/sigma)**3 )
    DZskew *= coef/DZstd**3
    Ristd = np.sqrt(coef*Ristd)    
    EAPE *= 0.5*coef

    return lon_deg, lat_deg, CTstd, SAstd, DZmean, DZstd, DZskew, Ristd, EAPE, NBstd


def main(itile, typestat, reso, timeflag):
    """Main function of stats.py"""
    create_stat_file(itile, typestat, reso, timeflag)
    write_stat_file(itile, typestat, reso, timeflag)
    #  create_stat_file(68, 'zmean', 0.5, 'annual')
    #  write_stat_file(68, 'zmean', 0.5, 'annual')


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    itile = 80
    reso_deg = 0.5
    timeflag = 'annual'
    lon_deg, lat_deg, CTstd, SAstd, DZmean, DZstd, DZskew, Ristd, EAPE, NBstd = compute_std_at_zref(itile, reso_deg, timeflag, verbose=True)
    #main(50, 'zstd', 0.5, 'annual')
#==============================================================================
#     main(51, 'zmean', 0.5, 'annual')
#     main(50, 'zmean', 0.5, 'annual')
#     main(70, 'zmean', 0.5, 'annual')
#     main(71, 'zmean', 0.5, 'annual')
#     main(72, 'zmean', 0.5, 'annual')
#==============================================================================
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
