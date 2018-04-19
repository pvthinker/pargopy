from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argodb as argo
import research_tools as research

plt.ion()
plt.close('all')

dirtopo = '/datawork/fsi2/mars/DATA/BATHY/ETOPO2'
topofile = 'etopo2.nc'

dirtopo = '/net/alpha/exports/sciences/data/BATHYMETRIE/BATHYMETRIE'
topofile = 'ETOPO2v2c_f4.nc'
dirtile = '/net/libra/local/tmp/1/herry/tiles'

itile = 50

figsize = (9, 7)

reso = 0.5

argodic = research.read_argo_filter(itile)
minlon, maxlon, minlat, maxlat = argodic['LONMIN_NO_M'], argodic['LONMAX_NO_M'], argodic['LATMIN_NO_M'], argodic['LATMAX_NO_M']

lon = np.arange(reso*np.floor(minlon/reso), reso*np.floor(maxlon/reso)+reso, reso)
lat = np.arange(reso*np.floor(minlat/reso), reso*np.floor(maxlat/reso)+reso, reso)

#lon_deg, lat_deg = define_grid(minlon, maxlon, minlat, maxlat, reso_deg)


with Dataset('%s/%s' % (dirtopo, topofile)) as nc:
    z = nc.variables['z'][:,:]


dl0 = 1/30.  # 1/30deg for etopo
lontopo = np.arange(-180, 180+dl0, dl0)
lattopo = np.arange(-90, 90+dl0, dl0)


def get_idx_of_box(lontopo, lattopo, cell):
    minlon, maxlon, minlat, maxlat = cell
    ilon = [i for i, x in enumerate(lontopo) if (x>=minlon) and (x<=maxlon)]
    jlon = [j for j, x in enumerate(lattopo) if (x>=minlat) and (x<=maxlat)]
    return ilon[0], ilon[-1], jlon[0], jlon[-1]

domain = [minlon, maxlon, minlat, maxlat]
i0, i1, j0, j1 = get_idx_of_box(lontopo, lattopo, domain)


def average_topo_on_box(depth, cell):
    """ average high resolution depth array on cell """
    i0, i1, j0, j1 = get_idx_of_box(lontopo, lattopo, cell)
    return np.mean(depth[j0:j1, i0:i1].ravel())


def box(cell, d=0):
    x1, x2, y1, y2 = cell
    plt.plot([x1-d, x1-d, x2+d, x2+d, x1-d],
           [y1-d, y2+d, y2+d, y1-d, y1-d], 'k')

plt.figure(figsize=figsize)
plt.imshow(z[j0:j1, i0:i1],
           origin='lower', extent=[minlon, maxlon, minlat, maxlat])
plt.axis('tight')
plt.colorbar()

reso = 0.5
lon = np.arange(minlon, maxlon, reso)
lat = np.arange(minlat, maxlat, reso)
    
nlon = len(lon)
nlat = len(lat)
bathy = np.zeros((nlat, nlon))
for j in range(nlat-1):
    for i in range(nlon-1):
        reso2 = reso*0.5
        gridcell = [lon[i]-reso2, lon[i]+reso2, lat[j]-reso2, lat[j]+reso2]
        box(gridcell)
        get_idx_of_box(lontopo, lattopo, gridcell)
        bathy[j, i] = average_topo_on_box(z, gridcell)

msk = bathy < 0

fig, ax = plt.subplots(2,1)

divider = make_axes_locatable(ax[0])
ax_cb = divider.new_horizontal(size="4%", pad=0.2)
im = ax[0].imshow(bathy,
           origin='lower', interpolation='nearest',
           extent=[minlon, maxlon, minlat, maxlat])
ax[0].set_title('tile #%03i' % itile)
fig.add_axes(ax_cb)
fig.colorbar(im, cax=ax_cb)

divider = make_axes_locatable(ax[1])
ax_cb = divider.new_horizontal(size="4%", pad=0.2)
ax[1].imshow(msk,
           origin='lower', interpolation='nearest',
           extent=[minlon, maxlon, minlat, maxlat])
