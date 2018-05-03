"""
Compute statistics on one tile

"""

from netCDF4 import Dataset
import os
import numpy as np
import gsw as gsw
import time
from scipy import interpolate as ip
import general_tools as tools
import argotools as argotools
import param as param
import melted_functions as melted
import netCDF_form as ncform

path_to_stats = param.path_to_stats
path_to_tiles = param.path_to_tiles
path_to_filter = param.path_to_filter
key_extraction = ['JULD', 'LONGITUDE', 'DATA_MODE', 'TAG', 'RHO', 'LATITUDE', 'SA', 'CT', 'BVF2']
zref = argotools.zref
daclist = argotools.daclist


def create_stat_file(itile, typestat, reso, timeflag, date, mode):
    """Create statistics netcdf file
    
    :rtype: None
    """
    
    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    argodic = melted.read_dic('argodic%003i' % itile, path_to_filter)
    minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']
    lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso)
    nlat, nlon = np.shape(lon_deg)

    ncform.netCDF_dim_creation(filename, zref, nlat, nlon, mode, date)
    ncform.netCDF_var_creation(filename, typestat)


def write_stat_file(itile, typestat, reso, timeflag, date, mode):
    """Write statistics into a netcdf file
    
    :rtype: None
    """

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    res = {}

    if typestat == 'zmean':
        res = compute_mean_at_zref(itile, reso, mode, date)
    elif typestat == 'zstd':
        res = compute_std_at_zref(itile, reso, timeflag, mode, date)

    res['zref'] = zref
    ncform.netCDF_var_writing(filename, typestat, res)


def read_stat_file(itile, typestat, reso, timeflag, date, mode):
    """Read statistics into a netcdf file
    
    :rtype: dict"""

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    print('read stat file : %s' % filename)

    res = ncform.netCDF_var_reading(filename, typestat)
    print(res)
    return res


def generate_filename(itile, typestat, reso, timeflag, date, mode):
    """
    Generates the filename of the netCDF stat file
    
    :rtype: str
    """
    filename = '%s/%s/%s/%s/%s/%s_%s_%s_%003i.nc' % (path_to_stats, 
                                                     reso, date[0], 
                                                     mode, typestat, 
                                                     typestat, reso, 
                                                     timeflag, itile)
    return filename


def grid_coordinate(itile, reso):
    """ 
    Returns the coordinates of each point of the grid for a given tile
    
    :rtype: numpy.ndarray, numpy.ndarray""" 
    lat, lon, nlat, nlon, marginlat, marginlon = argotools.tile_definition()

    i = itile / nlon
    j = itile - (i * nlon)

    latmin = lat[i]
    latmax = lat[i+1]
    lonmin = lon[j]
    lonmax = lon[j+1]

    grid_lat = np.arange(latmin, (latmax+reso), reso)
    grid_lon = np.arange(lonmin, (lonmax+reso), reso)

    return(grid_lat, grid_lon)
    


def define_grid(minlon, maxlon, minlat, maxlat, reso_deg):
    """ setup the grid coordinates (in degrees)
    coordinates are round multiples of reso_deg
    reso_deg sets the grid resolution, typically 0.5deg
    
    :rtype: float, float""" 

    minlon = np.ceil(minlon/reso_deg)*reso_deg
    maxlon = np.floor(maxlon/reso_deg)*reso_deg

    minlat = np.ceil(minlat/reso_deg)*reso_deg
    maxlat = np.floor(maxlat/reso_deg)*reso_deg

    lon1D = np.arange(minlon, maxlon+reso_deg, reso_deg)
    lat1D = np.arange(minlat, maxlat+reso_deg, reso_deg)

    lon, lat = np.meshgrid(lon1D, lat1D)

    return lon, lat


def compute_mean_at_zref(itile, reso_deg, mode, date):
    """Compute the mean at depths zref
    
    :rtype: dict"""
    tile = data_choice(mode, date, itile)
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']
    #  argodic = argotools.read_argo_filter(itile)
    argodic = melted.read_dic('argodic%003i' % itile, path_to_filter)
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
        #  print('%4i/%i' % (k, nbprof))
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

    res = {'NBbar': NBbar,
           'CTbar': CTbar,
           'SAbar': SAbar,
           'Ribar': RIbar,
           'BVF2bar': BVF2bar,
           'lon': lon_deg,
           'lat': lat_deg}

    return res


def compute_std_at_zref(itile, reso_deg, timeflag, mode, date, verbose=False):
    """Compute the standard deviations at depths zref
    
    :rtype: dict"""

    # gridded arrays of CT, SA variances
    res = read_stat_file(itile, 'zmean', reso_deg, timeflag, date, mode) # read it from the file
    tile = data_choice(mode, date, itile)
    # output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    # iprof = output['IPROF']
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']
    
    lon_rad = np.deg2rad(res['lon'])
    lat_rad = np.deg2rad(res['lat'])
    reso_rad = np.deg2rad(reso_deg)

    nlat, nlon = np.shape(res['lon'])

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
            interpolator = ip.interp1d(res['Ribar'][:,j,i], zref, bounds_error=False)
            p = gsw.p_from_z(-zref, lat[j])
            g = gsw.grav(lat[j], p)
            cs = gsw.sound_speed(res['SAbar'][:, j, i], res['CTbar'][:, j, i], p)
            rho0 = res['Ribar'][:, j, i].copy()

            drho = RI - res['Ribar'][:, j, i]
            dbvf2 = BVF2 - res['BVF2bar'][:, j, i]
            dCT = CT - res['CTbar'][:, j, i]
            dSA = SA - res['SAbar'][:, j, i]
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
            BVF2std[:, j, i] = average(dbvf2**2)
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
    BVF2std = np.sqrt(coef*BVF2std)
    EAPE *= 0.5*coef

    res = {'NBstd': NBstd,
           'CTstd': CTstd,
           'SAstd': SAstd,
           'Ristd': Ristd,
           'BVF2std': BVF2std,
           'DZmean': DZmean,
           'DZstd': DZstd,
           'DZskew': DZskew,
           'EAPE': EAPE,
           'lon': res['lon'],
           'lat': res['lat']}

    return res


def data_choice(mode, date, itile):
    """
    Make the tile filter to choose keep only the chosen mode ('R', 'A', 'D', 'AD' or 'RAD')
    and the profiles under the chosen date (year, month, day)
    Return the tile_extract according to the filters used.
    
    :rtype: dic
    """
    tile = melted.read_dic('tile%003i' % itile, path_to_tiles)
    #  tile = tiler.read_tile(itile)
    julday = argotools.conversion_gregd_juld(int(date[0]), int(date[1]), int(date[2]))
    mode_list = list(mode)
    idx = []
    tile_extract = {}
    for m in mode_list:
            idx += np.where((tile['DATA_MODE'] == m) & (tile['JULD'] < julday))
    if len(mode_list) == 2:
        idx1 = np.concatenate((idx[0], idx[1]))
    elif len(mode_list) == 3:
        idx1 = np.concatenate((idx[0], idx[1], idx[2]))
    else:
        idx1 = idx[0]
    for k in key_extraction:
        tile_extract[k] = tile[k][idx1]
    tile_extract['ZREF'] = tile['ZREF']
        

    return tile_extract


def main(itile, typestat, reso, timeflag, date, mode):
    """Main function of stats.py"""
    create_stat_file(itile, typestat, reso, timeflag, date, mode)
    write_stat_file(itile, typestat, reso, timeflag, date, mode)
    #  grid_coordinate(299, 0.5)
    #  read_stat_file(itile, typestat, reso, timeflag, date, mode)
    #  create_stat_file(68, 'zmean', 0.5, 'annual')
    #  write_stat_file(68, 'zmean', 0.5, 'annual')


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    main(32, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'D')
    #  main(293, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
    #  main(294, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
    #  main(295, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
    #  main(296, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
    #  main(297, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
    #  main(298, 'zstd', 0.5, 'annual', ['2017', '12', '31'], 'AD')
#==============================================================================
#     main(51, 'zmean', 0.5, 'annual')
#     main(50, 'zmean', 0.5, 'annual')
#     main(70, 'zmean', 0.5, 'annual')
#     main(71, 'zmean', 0.5, 'annual')
#     main(72, 'zmean', 0.5, 'annual')
#==============================================================================
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)