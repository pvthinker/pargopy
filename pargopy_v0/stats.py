"""
Compute statistics on one tile

"""

import os
import numpy as np
import gsw as gsw
import time
from scipy import interpolate as ip
import general_tools as general
import argotools as argotools
import param as param
import netCDF_form as ncform
import variable_selector as vs
import itertools as it
import atlas as atlas
import eape as eape

# Put in var_choice list the variables you want in your atlas

#  var_choice = ['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std', 'DZmean', 'DZstd', 'DZskew', 'EAPE']
block_choice = ['zmean', 'zstd', 'zdz']  # , 'rmean', 'rstd', 'rdz'
var_choice = {'zmean': ['NBbar', 'SAbar', 'CTbar', 'Ribar', 'BVF2bar'],
              'zstd': ['NBstd', 'SAstd', 'CTstd', 'Ristd', 'BVF2std'], 
              'zdz': ['DZmean', 'DZstd', 'DZskew', 'EAPE']}
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

    print('create stat file %s' % filename)
    ncform.create_dim(filename, zref, nlat, nlon, mode, date)
    ncform.create_var(filename, var_choice[typestat])


def write_stat_file(itile, typestat, reso, timeflag, date, mode, stats_mode):
    """Write statistics into a netcdf file
    
    :rtype: None
    """

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    res = {}

    res, nanidx = compute_at_zref(itile, reso, mode, date, stats_mode)

    res['zref'] = zref
    ncform.write_var(filename, var_choice[typestat], res)

    filename_zmean = generate_filename(itile, 'zmean', reso, timeflag, date, mode)
    if os.path.isfile(filename_zmean):
        print('File of zmean stats already written')
    else:
        print('File of zmean stats not write, let\'s do it...')
        ncform.write_var(filename_zmean, var_choice['zmean'], res)
        print('File of zmean stats written')


def read_stat_file(itile, typestat, reso, timeflag, date, mode, var_choice):
    """Read statistics into a netcdf file
    
    :rtype: dict"""

    filename = generate_filename(itile, typestat, reso, timeflag, date, mode)
    print('read stat file : %s' % filename)

    res = ncform.read_var(filename, var_choice)
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
    


def compute_at_zref(itile, reso_deg, mode, date, block_choice, tile_dict=None):
    """Compute the variables at depths zref
    
    :rtype: dict"""
    if tile_dict != None:
        tile= tile_dict
    else:
        tile = date_mode_filter(mode, date, itile)
    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    nanidx = np.where(np.isnan(CT) | np.isnan(SA))
    print(nanidx)
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

    variables = {}

    for b in block_choice:
        # gridded arrays of CT, SA et RI means
        for i, v in enumerate(var_choice['zmean']):
            variables[v] = np.zeros((nz, nlat, nlon))

        for k in range(nbprof):
            #  print('%4i/%i' % (k, nbprof))
            # todo: weigh in time using juld,
            # e.g. only winter statistics
            time_weight = 1.
            xlon_rad = np.deg2rad(lon[k])
            xlat_rad = np.deg2rad(lat[k])
            weight = general.compute_weight(lon_rad, lat_rad,
                                            xlon_rad, xlat_rad,
                                            reso_rad)
            weight *= time_weight
            for l in range(nz):
                if np.isnan(CT[k, l]) or np.isnan(SA[k, l]):
                    pass
                else:
                    variables['NBbar'][l, :, :] += weight
                    variables['CTbar'][l, :, :] += weight*CT[k, l]
                    variables['SAbar'][l, :, :] += weight*SA[k, l]
                    variables['Ribar'][l, :, :] += weight*RI[k, l]
                    variables['BVF2bar'][l, :, :] += weight*BVF2[k, l]

        # normalize with the number of profiles (fractional
        # because NBbar is fractionnal)
        coef = 1./variables['NBbar']
        coef[variables['NBbar'] < 1] = np.NaN

        variables['CTbar'] *= coef
        variables['SAbar'] *= coef
        variables['Ribar'] *= coef
        variables['BVF2bar'] *= coef


        if b =='zstd' or b == 'zdz':
            xlon_rad = np.deg2rad(lon)
            xlat_rad = np.deg2rad(lat)
            for i, v in enumerate(var_choice[b]):
                variables[v] = np.zeros((nz, nlat, nlon))
            variables['NBstd'] = variables['NBbar']

            if len(lat) == 0:
                pass
            else:
                for j in range(nlat):
                    for i in range(nlon):
                        if len(lat) < j+1:
                            pass
                        else:
                            time_weight = 1.
                            weight = general.compute_weight(lon_rad[j, i], lat_rad[j, i],
                                                            xlon_rad, xlat_rad,
                                                            reso_rad)
                            weight *= time_weight
                            drho = RI - variables['Ribar'][:, j, i]
                            dbvf2 = BVF2 - variables['BVF2bar'][:, j, i]
                            dCT = CT - variables['CTbar'][:, j, i]
                            interpolator = ip.interp1d(variables['Ribar'][:,j,i], zref, bounds_error=False)
                            p = gsw.p_from_z(-zref, lat[j])
                            g = gsw.grav(lat[j], p)
                            cs = gsw.sound_speed(variables['SAbar'][:, j, i], variables['CTbar'][:, j, i], p)
                            rho0 = variables['Ribar'][:, j, i].copy()
                            zrho = interpolator(RI)
                            dzstar = zrho-zref
                            dz = dzstar/(1.+rho0*g*dzstar/(cs**2*drho))
                            dSA = SA - variables['SAbar'][:, j, i]

                            weight = weight[:, np.newaxis] + np.zeros_like(zref)
                            weight[np.where(np.isnan(dz) | np.isnan(drho) | np.isnan(dCT) | np.isnan(dSA))] = 0.
                            weight[nanidx] = 0.
                            def average(field):
                                return np.nansum(weight*field, axis=0)
                            if b == 'zstd':
                                variables['CTstd'][:, j, i] = average(dCT**2)
                                variables['SAstd'][:, j, i] = average(dSA**2)
                                variables['BVF2std'][:, j, i] = average(dbvf2**2)
                                variables['Ristd'][:, j, i] = average(drho**2)
                               
                                coef = 1./(variables['NBstd']-1)
                                coef[variables['NBstd'] < 2] = np.nan

                                

                            if b == 'zdz':
            
                                variables['DZmean'][:, j, i] = average(dz)
                                variables['DZstd'][:, j, i] = average(dz**2)
                                variables['DZskew'][:, j, i] = average(dz**3)
                                variables['EAPE'][:, j, i] = average(dz*drho)
                                
                                
        if b == 'zstd':
            variables['CTstd'] = np.sqrt(coef*variables['CTstd'])
            variables['SAstd'] = np.sqrt(coef*variables['SAstd'])
            variables['Ristd'] = np.sqrt(coef*variables['Ristd'])
            variables['BVF2std'] = np.sqrt(coef*variables['BVF2std'])
            
        elif b == 'zdz':
            variables['DZmean'] *= coef
            variables['DZstd'] = np.sqrt(coef*variables['DZstd'])
            variables['DZskew'] *= coef/variables['DZstd']**3
            variables['EAPE'] *= 0.5*coef


    variables['lat'] = lat_deg
    variables['lon'] = lon_deg
    return variables, nanidx, weight

def retrieve_tile_from_position(lon0, lat0):
    """Return the tile index in which (lon0, lat0) sits
    
    :rtype: list"""

    lat, lon, nlat, nlon, marginlat, marginlon = argotools.tile_definition()
    j = [k for k in range(len(lat)) if lat[k]<lat0][-1]
    i = [k for k in range(len(lon)) if lon[k]<lon0][-1]
    return atlas.ij2tile(i, j)


def compute_stats_at_zref(mode, date, grid_lon, grid_lat, reso_deg, manual_check=False, new_tile=None):
    """ compute statistics on a small grid defined at grid_lon x grid_lat
    the small grid should fit inside one tile """
    tmps1 = time.time()
    itiles = [retrieve_tile_from_position(lon0, lat0) for lon0,lat0 in it.product(grid_lon,grid_lat)]
    if len(set(itiles)) == 1:
        itile = itiles[0]
        if manual_check == True:
            tile = new_tile
        else:
            tile = date_mode_filter(mode, date, itile)
    else:
        raise ValueError('the domain does not fit in one tile')

    CT, SA, RI, BVF2 = tile['CT'], tile['SA'], tile['RHO'], tile['BVF2']
    lat, lon = tile['LATITUDE'], tile['LONGITUDE']

    xlon_rad = np.deg2rad(lon)
    xlat_rad = np.deg2rad(lat)
    reso_rad = np.deg2rad(reso_deg)

    nlon = len(grid_lon)
    nlat = len(grid_lat)
    nz = len(zref)

    Nb = np.zeros((nz, nlat, nlon))
    CTbar = np.zeros((nz, nlat, nlon))
    CTstd = np.zeros((nz, nlat, nlon))
    SAbar = np.zeros((nz, nlat, nlon))
    SAstd = np.zeros((nz, nlat, nlon))
    RIbar = np.zeros((nz, nlat, nlon))
    RIstd = np.zeros((nz, nlat, nlon))

    EAPE = np.zeros((nz, nlat, nlon))
    DZstd = np.zeros((nz, nlat, nlon))

    def average(field):
        """ weighted average of profiles

        field is a nprofiles x nz array
        the weight depends on the profile to grid point distance"""

        return np.nansum(weight*field, axis=0)

    def bar_std(phi, coef, coef1):
        """ compute the mean and the std of 'phi', where phi
        are profiles interpolated at zref of either CT, SA, RI etc."""

        phibar = coef*average(phi)
        dphi = phi - phibar
        phistd = np.sqrt(coef1*average(dphi**2))
        return phibar, phistd

    nanidx = np.where(np.isnan(CT) | np.isnan(SA))

    for j, lat_rad in enumerate(np.deg2rad(grid_lat)):
        for i, lon_rad in enumerate(np.deg2rad(grid_lon)):
            #print('%i/%i - %i/%i' % (i, nlon, j, nlat))
            weight = general.compute_weight(lon_rad, lat_rad,
                                          xlon_rad, xlat_rad, reso_rad)
            weight = weight[:, np.newaxis] + np.zeros_like(zref)
            weight[nanidx] = 0.

            Nb[:, j, i] = average(1.)
            coef = 1./Nb[:, j, i]
            coef[Nb[:, j, i]<2] = np.NaN
            coef1 = 1./(Nb[:, j, i]-1)
            coef1[Nb[:, j, i]<2] = np.NaN

            CTbar[:, j, i], CTstd[:, j, i] = bar_std(CT, coef, coef1)
            SAbar[:, j, i], SAstd[:, j, i] = bar_std(SA, coef, coef1)
            RIbar[:, j, i], RIstd[:, j, i] = bar_std(RI, coef, coef1)

            # the original piece of code, before defining bar_std()
            #
            # CTbar[:, j, i] = coef*average(CT)            
            # dCT = CT - CTbar[:, j, i]
            # CTstd[:, j, i] = np.sqrt(coef1*average(dCT**2))

            # RIbar[:, j, i] = coef*average(RI)            
            # dRI = RI - RIbar[:, j, i]
            # RIstd[:, j, i] = np.sqrt(coef1*average(dRI**2))

            # SAbar[:, j, i] = coef*average(SA)            
            # dSA = SA - SAbar[:, j, i]
            # SAstd[:, j, i] = np.sqrt(coef1*average(dSA**2))

            p = gsw.p_from_z(-zref, grid_lat[j])
            
            cs = gsw.sound_speed(SAbar[:, j, i], CTbar[:, j, i], p)

            rho0 = RIbar[:, j, i].copy()

            zr, Eape = eape.compute_eape(zref, rho0, cs, RI)

            dZ = zr-zref
            DZstd[:, j, i] = np.sqrt(coef1*average(dZ**2))
            EAPE[:, j, i] = np.sqrt(coef1*average(Eape))

    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)

    return {'NB': Nb,
            'lon': grid_lon, 'lat': grid_lat, 'zref': zref,
            'CTbar': CTbar, 'CTstd': CTstd,
            'RIbar': RIbar, 'RIstd': RIstd,
            'SAbar': SAbar, 'SAstd': SAstd,
            'DZstd': DZstd, 'EAPE': EAPE}



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
    main(52, ['zmean'], 0.5, 'annual', ['2017', '12', '31'], 'D')
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
