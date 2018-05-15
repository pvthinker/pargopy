"""
Compute statistics on one tile

"""

import os
import numpy as np
import gsw as gsw
import time
from scipy import interpolate as ip
import general_tools as tools
import argotools as argotools
import param as param
import netCDF_form as ncform
import variable_selector as vs

# Put in var_choice list the variables you want in your atlas

#  var_choice = ['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std', 'DZmean', 'DZstd', 'DZskew', 'EAPE']
var_choice = {'zmean':['NBbar', 'SAbar', 'CTbar', 'Ribar', 'BVF2bar'],
              'zstd':['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std'], 
              'zdz':['DZmean', 'DZstd', 'DZskew', 'EAPE']}
# var_choice used by stats_read
# You need to know which variables are existing in the stats_i.nc file to know
# which variables are available
# !-- to do : Replace the reading of the stats_i.nc file by calculating the
#             XXbar variables needed to calculate the XXstd variables wanted --!
var_choice_mean = ['NBbar', 'SAbar', 'CTbar', 'Ribar', 'BVF2bar']

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
    grid_lat, grid_lon = grid_coordinate(itile, reso)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)
    nlat, nlon = np.shape(lon_deg)

    print(filename)
    ncform.netCDF_dim_creation(filename, zref, nlat, nlon, mode, date)
    ncform.netCDF_var_creation(filename, var_choice[typestat])


def write_stat_file(itile, typestat, reso, timeflag, date, mode, stats_mode):
    """Write statistics into a netcdf file
    
    :rtype: None
    """

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    res = {}

    res = vs.compute_at_zref(itile, reso, mode, date, stats_mode)

#==============================================================================
#     if typestat[0] == 'zmean':
#         res = compute_mean_at_zref(itile, reso, mode, date)
#     elif typestat[0] == 'zstd':
#         res = compute_std_at_zref(itile, reso, timeflag, mode, date)
#==============================================================================

    res['zref'] = zref
    ncform.netCDF_var_writing(filename, var_choice[typestat], res)

    filename_zmean = generate_filename(itile, 'zmean', reso, timeflag, date, mode)
    if os.path.isfile(filename_zmean):
        print('File of zmean stats already written')
    else:
        print('File of zmean stats not write, let\'s do it...')
        ncform.netCDF_var_writing(filename_zmean, var_choice['zmean'], res)
        print('File of zmean stats written')


def read_stat_file(itile, typestat, reso, timeflag, date, mode, var_choice):
    """Read statistics into a netcdf file
    
    :rtype: dict"""

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    print('read stat file : %s' % filename)

    res = ncform.netCDF_var_reading(filename, var_choice)
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

    coordinates are round multiples of reso_deg
    reso sets the grid resolution, typically 0.5deg
    
    :rtype: numpy.ndarray, numpy.ndarray""" 
    lat, lon, nlat, nlon, marginlat, marginlon = argotools.tile_definition()

    i = itile // nlon
    j = itile % nlon

    latmin = np.ceil(lat[i]/reso)*reso
    latmax = np.floor(lat[i+1]/reso)*reso
    lonmin = np.ceil(lon[j]/reso)*reso
    lonmax = np.floor(lon[j+1]/reso)*reso
    if (lonmax % 1. == 0.):
        # remove rightmost point (it's the leftmost point of the next
        # tile). That's also true for +180 (that is == -180).
        lonmax -= reso

    grid_lat = np.arange(latmin, (latmax+reso), reso)
    grid_lon = np.arange(lonmin, (lonmax+reso), reso)

    return(grid_lat, grid_lon)
    


def compute_mean_at_zref(itile, reso_deg, mode, date):
    """Compute the mean at depths zref
    
    :rtype: dict"""
    tile = date_mode_filter(mode, date, itile)
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']
    grid_lat, grid_lon = grid_coordinate(itile, reso_deg)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)

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

    #  print(CTbar)
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
    res = read_stat_file(itile, 'zmean', reso_deg, timeflag, date, mode, var_choice_mean) # read it from the file
    tile = date_mode_filter(mode, date, itile)
    # output = argotools.retrieve_infos_from_tag(argodb, tile['TAG'])
    # iprof = output['IPROF']
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']

    grid_lat, grid_lon = grid_coordinate(itile, reso_deg)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)
    nlat, nlon = np.shape(lon_deg)

    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)
    reso_rad = np.deg2rad(reso_deg)

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
    if len(lat) == 0:
        pass
    else:
        for j in range(nlat):
            for i in range(nlon):
                if len(lat) < j+1:
                    pass
                else:
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
                    
                    # zstd
                    CTstd[:, j, i] = average(dCT**2)
                    SAstd[:, j, i] = average(dSA**2)
                    BVF2std[:, j, i] = average(dbvf2**2)
                    Ristd[:, j, i] = average(drho**2)

                    # zdz
                    DZmean[:, j, i] = average(dz)
                    DZstd[:, j, i] = average(dz**2)
                    DZskew[:, j, i] = average(dz**3)
                    EAPE[:, j, i] = average(dz*drho)

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


def date_mode_filter(mode, date, itile):
    """
    Make the tile filter to choose keep only the chosen mode ('R', 'A', 'D', 'AD' or 'RAD')
    and the profiles under the chosen date (year, month, day)
    Return the tile_extract according to the filters used.
    
    :rtype: dic
    """
    tile = argotools.read_dic('tile%003i' % itile, path_to_tiles)
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
    for t in typestat:
        create_stat_file(itile, t, reso, timeflag, date, mode)
        write_stat_file(itile, t, reso, timeflag, date, mode, typestat)
    #  grid_coordinate(itile, reso)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    main(52, ['zdz', 'zstd'], 0.5, 'annual', ['2017', '12', '31'], 'D')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)