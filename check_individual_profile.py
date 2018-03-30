from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from tile import Tile
import general_tools as gr
import gsw as gsw
import profile_tools as pr
import research_tools as res
import interpolation_tools as interpolation

#plt.ion()

dirtile = '/net/libra/local/tmp/1/herry/tiles'
dirlog = '/net/libra/local/tmp/1/herry/data'

# nc = Dataset('%s/research.nc' % (dirlog))
# nc.close()


itile = 51

nc = Dataset('%s/tile%03i.nc' % (dirtile, itile))
nb_prof = len(nc.dimensions['nb_prof'])

minlon = nc.lonmin
maxlon = nc.lonmax
minlat = nc.latmin
maxlat = nc.latmax
print('lon :', minlon, maxlon)
print('lat :', minlat, maxlat)

CT = nc.variables['CT'][:, :]
SA = nc.variables['SA'][:, :]
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
prof = nc.variables['prof'][:]
juld = nc.variables['juld'][:]
wmo = nc.variables['wmo'][:]
dac = nc.variables['dac'][:]
zref = nc.variables['zref'][:]
nc.close()
nbprof = len(prof)
nz = len(zref)

# profile index
ok = False
l = 100
profok = []
ctl = np.sum(np.isnan(CT), axis=1)
for l in range(nbprof):
    if ctl[l] == 0:
        profok.append(l)
l = np.random.choice(profok, size=1)[0]
#l = 981
print('found profile %i' % l)

temp, sal, pres, temp_qc, sal_qc, pres_qc, lat0, lon0, juld0 = pr.read_profile(dac[l], wmo[l])


plt.figure()
n, nz = np.shape(sal)
for idx in range(3):
    Ti, Si, Ri, ct, sa, z, ierr = interpolation.raw_to_interpolate(
        temp[idx], sal[idx], pres[idx],
        temp_qc[idx], sal_qc[idx], pres_qc[idx],
        lon0[idx], lat0[idx], zref)

    if ierr==0:
        print(len(sa), len(z))
        plt.plot(sa, -z, 'k+-')
        plt.plot(Si, -zref, 'bo-', lw=2)

# plt.plot(SA[l, :], -zref,'r')
plt.show()
iprof=prof[l]
print(lon0[iprof],lat0[iprof],juld0[iprof])
print(lon[l],lat[l],juld[l])
