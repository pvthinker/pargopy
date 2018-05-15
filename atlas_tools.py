"""
Created on Mon Mar 12 13:10:24 2018

File used to create the on-click tools for the atlas

"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4  import Dataset
import research_tools as res
import atlas as at
import general_tools as gene
import param as param
import argotools as argotools
import mouseprofile as mouse
import tempfile

diratlas = param.path_to_atlas
dirtile = param.path_to_tiles

atlas_name = at.atlas_name
listvar = at.listvar

reso = 0.5
nlon = 20
nlat = 15

#plt.ion()

def lonlatstr(lon, lat):
    if lon > 0:
        lonchar = 'E'
    else:
        lonchar = 'W'
    if lat > 0:
        latchar = 'N'
    else:
        latchar = 'S'
    pos = '%.2f%s - %.2f%s' % (abs(lon), lonchar, abs(lat), latchar)
    return pos

def retrieve_tile_from_position(lon0, lat0):
    """Return the tile index in which (lon0, lat0) sits
    
    :rtype: list"""

    lat, lon, nlat, nlon, marginlat, marginlon = argotools.tile_definition()
    j = [k for k in range(len(lat)) if lat[k]<lat0][-1]
    i = [k for k in range(len(lon)) if lon[k]<lon0][-1]
    return at.ij2tile(i, j)


def select_profiles_near_point(lon0, lat0):
    """Return the list of profiles at a distance 'reso' from
    (lon0, lat0)
    
    :rtype: None"""
    
    itile = retrieve_tile_from_position(lon0, lat0)
    argodb = argotools.read_dic('tile%003i' % itile, dirtile)
    zref = argodb.pop('ZREF')
    dmin = np.deg2rad(reso)
    lonpro = argodb['LONGITUDE']
    latpro = argodb['LATITUDE']
    dist = gene.dist_sphe(np.deg2rad(lonpro), np.deg2rad(latpro),
                          np.deg2rad(lon0), np.deg2rad(lat0))
    idx = np.where(dist < dmin)[0]
    
    subargo = res.extract_idx_from_argodb(argodb, idx)

    #fig2 = plt.figure(2 ,figsize=(16,6))
    fig2.clf()
    var = ['CT', 'SA', 'RHO']
    var = ['SA']
    for j, v in enumerate(var):
        ax = fig2.add_subplot(111+j)
        if v == 'SA':
             for k in range(len(mps)):
                 # empty this list of graphical objects
                 # to prevent it from getting too fat
                 mps.pop()
                 
        for k in range(len(subargo['JULD'])):
            li, = ax.plot(subargo[v][k, :], -zref, picker=5)
            ax.set_xlabel(v)
            ax.set_ylabel('z [m]')
            ax.set_title('profiles near %s / tile = %i' %
                      (lonlatstr(lon0, lat0), itile))
            if v == 'SA':
                tag = subargo['TAG'][k]
                mp = mouse.MouseProfile(li, tag, badprofiles)
                mp.connect()
                mps.append(mp)
    fig2.canvas.draw()
    return

# ncfile = '%s/%s.nc' % (diratlas, atlas_name)
ncfile = '/home2/datawork/therry/tmp/atlas/0.5/2017/D/zstd/zstd_0.5_annual.nc'

with Dataset(ncfile, 'r', format='NETCDF4') as nc:
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    # zref = nc.variables['zref'][:]
    EAPE =  nc.variables['EAPE'][:, :, :]

longrid, latgrid = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))
#CT.mask = False

# plt.ion()

zref = at.zref
#fig2 = plt.figure(2)
fig2 = plt.figure(2 ,figsize=(16,6))


# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print('Before gridindex')
    lon, lat = at.gridindex2lonlat(ix, iy)
    print('After gridindex')
    itile = retrieve_tile_from_position(lon, lat)
    print('x = %d, y = %d, lon = %d, lat = %d, tile = %i'
          % (ix, iy, lon, lat, itile))


    select_profiles_near_point(lon, lat)
    
    # assign global variable to access outside of function

    coords.append((ix, iy))

    # Disconnect after 2 clicks
    if len(coords) == 2:
        pass
        #fig.canvas.mpl_disconnect(cid)
        #plt.close(1)
    return

global coords
global mps
mps = []
coords = []

spoiled_profiles=[]

badprofiles = 'badprofiles.txt'

fig = plt.figure(1)
kz = 20

# plt.pcolor(lon, lat, CT[kz, :, :])
plt.imshow(EAPE[kz, :,:], origin='lower',
           vmin=0, vmax=10.,
           interpolation='nearest', 
           cmap=plt.get_cmap('RdBu_r'))
plt.title('EAPE @ z=%.0fm' % zref[kz])
plt.colorbar()
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()