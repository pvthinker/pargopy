import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.interpolate import interp1d
from scipy.integrate import fixed_quad
from interp import PLagrange
import tile as tile
import gsw as gsw
from general_tools import dist_sphe
import warnings


warnings.filterwarnings("ignore", category=RuntimeWarning)


def wherenotnan(x):
    return [k for k, xx in enumerate(x) if not(np.isnan(xx))]


def wherenotnanxy(x, y):
    return [k for k in range(len(x)) if not(np.isnan(x[k])) and not(np.isnan(y[k]))]


def compute_eape(z0, r0, cs, rho):
    """Compute the eape of rho(...,z0) based on reference profile r0(z0) with
    sound_speed cs(z0)"""
    g = 9.81

    zr = np.zeros_like(rho)
    eape = np.zeros_like(rho)

    rhorefinterp = PLagrange(z0, r0, order=4)
    primitive = rhorefinterp.antiderivative()

    for k in range(len(z0)):
        cff = g*r0[k]/cs[k]**2

        def h(x): return rhorefinterp(x) - cff*(x-z0[k])

        zzr = interp1d(h(z0), z0, bounds_error=False, copy=True,
                       kind='linear', assume_sorted=True)

        def area(x): return primitive(x)-cff*(0.5*x-z0[k])*x

        zr[..., k] = zzr(rho[..., k])
        eape[..., k] = (-area(zr[..., k])+area(z0[k])
                        - rho[..., k]*(z0[k]-zr[..., k]))

    return zr, eape


def integrate_K(zref, z, K):
    """ compute \int_0^z K(z')dz'"""

    zm = 0.5*(zref[1:]+zref[:-1])
    dz = zref[1:]-zref[:-1]
    Kref = interp1d(z, K, bounds_error=False, kind='linear',
                    fill_value='extrapolate')(zm)

    intKref = zref*0.
    intKref[1:] = np.cumsum(Kref*dz)
    intK = [intKref[k] for k, zz in enumerate(zref) if zz in z]

    return np.asarray(intK)

if __name__ == '__main__':

    plt.ion()

    tile.path_to_tiles = "/home/roullet/data/ARGO/tiles"
    zrefprofiles = tile.read_tile(51)

    nprof = len(zrefprofiles[b'JULD'])
    i0 = np.random.randint(nprof)

    # i0, i1 = 12216, 6482
    # i0, i1 = 8859, 10559  # pb at 700m out of the blue

    # i0, i1 = (11427, 3576)
    # i0, i1 = (2004, 2000)

    zref = zrefprofiles[b'ZREF']
    rhoref = zrefprofiles[b'RHO'][i0, :]
    CT = zrefprofiles[b'CT'][i0, :]
    SA = zrefprofiles[b'SA'][i0, :]
    lat = zrefprofiles[b'LATITUDE'][i0]
    lon = zrefprofiles[b'LONGITUDE'][i0]
    p0 = gsw.p_from_z(-zref, lat)
    # g = gsw.grav(lat, p)
    # cs = gsw.sound_speed(SA, CT, p0)

    lat1 = zrefprofiles[b'LATITUDE'][:]
    lon1 = zrefprofiles[b'LONGITUDE'][:]

    dist = dist_sphe(np.deg2rad(lon), np.deg2rad(lat),
                     np.deg2rad(lon1), np.deg2rad(lat1))

    dist = np.rad2deg(dist)

    dmax = 1.

    idx = list(np.where(dist < dmax)[0])
    idxpos = idx.copy()
    # idx.remove(i0)
    CT = np.nanmean(zrefprofiles[b'CT'][idx, :], axis=0)
    SA = np.nanmean(zrefprofiles[b'SA'][idx, :], axis=0)
    rhoref = gsw.rho(SA, CT, p0)
    cs = gsw.sound_speed(SA, CT, p0)

    i1 = np.random.choice(idx)


    rho = zrefprofiles[b'RHO'][i1, :]
    CT = zrefprofiles[b'CT'][i1, :]
    SA = zrefprofiles[b'SA'][i1, :]
    lat1 = zrefprofiles[b'LATITUDE'][i1]
    lon1 = zrefprofiles[b'LONGITUDE'][i1]
    p = gsw.p_from_z(-zref, lat1)
    # g = gsw.grav(lat, p)
    cs1 = gsw.sound_speed(SA, CT, p)
    # print(rhoref, rho)

    # plt.plot(rhoref, -zref)
    # plt.plot(rho, -zref)

    dist = dist_sphe(np.deg2rad(lon), np.deg2rad(lat),
                     np.deg2rad(lon1), np.deg2rad(lat1))

    dist = np.rad2deg(dist)
    print('dist = %.2f' % dist)
    idx = wherenotnanxy(rhoref, rho)

    if idx and (dist < dmax):
        zr3, eape3 = compute_eape(zref[idx], rhoref[idx], cs[idx], rho[idx])

        K = 1/cs[idx]**2  # 2e-7  # should be 1/cs^2
        K1 = 1/cs1[idx]**2
        Km = 0.5*(K+K1)
        # g = g[0]
        g = 9.81

        z0 = zref[idx]
        r0 = rhoref[idx]
        rho = rho[idx]

        # rhor = interp1d(z0, r0, bounds_error=False,
        #                 kind='cubic', fill_value='extrapolate')
        rhor = PLagrange(z0, r0, order=4)
        # zzr0 = PLagrange(r0, z0, order=4)
        # determine zr
        zr = np.zeros_like(rho)
        zr2 = np.zeros_like(rho)
        eape = np.zeros_like(rho)
        eape2 = np.zeros_like(rho)
        plt.figure(1)
        plt.clf()
        # plt.plot(rhor(z0)-1020., -z0, lw=2)

        r0star = r0*np.exp(-g*z0*K)
        rhostar = rho*np.exp(-g*z0*K1)

        Kint = integrate_K(zref, z0, K)
        Kint1 = integrate_K(zref, z0, K1)

        r0star = r0*np.exp(-g*Kint)
        # il faut absolument prendre Kint et pas Kint1 !!!
        rhostar = rho*np.exp(-g*Kint)

        zstar = interp1d(r0star, z0, bounds_error=False,
                         kind='cubic')  # , fill_value='extrapolate')

        # zstar = PLagrange(r0star, z0, order=2)

        for k in range(len(z0)):
            def compress(x): return g*K[k]*r0[k]*(x-z0[k])

            def f(x): return rhor(x)-rho[k] - compress(x)

            def ff(x): return rhor(x)-rho[k]

            def h(x): return rhor(x) - compress(x)

            zzr = interp1d(h(z0), z0, bounds_error=False,
                           kind='cubic')

            # zzr = PLagrange(r0, z0, order=4)
            # def hh(x): return +(x-r0[k])/(g*K[k]*rho[k])

            # def zzr(x): return zzr0(x) - hh(x)
            try:
                # zr[k] = newton(f, z0[k], tol=1e-1, maxiter=7)
                zr[k] = +zstar(rhostar[k])
                zr2[k] = zzr(rho[k])

                eape[k], dummy = fixed_quad(f, zr[k], z0[k], n=8)
                eape2[k], dummy = fixed_quad(h, zr[k], z0[k], n=8)
                eape2[k] -= rho[k]*(z0[k]-zr[k])
                # eape2[k] = eape[k]
            except RuntimeError:
                zr[k] = np.nan
                zr2[k] = np.nan
                eape[k] = np.nan
                eape2[k] = np.nan

            if z0[k] == 1500:
                print(z0[k])
                plt.figure(3)
                plt.clf()
                plt.plot(z0, f(z0))
                plt.plot(z0, ff(z0))
                plt.plot(z0, z0*0, 'k')
                plt.plot([z0[k], z0[k]], [-1, 1], 'k')
                plt.plot(zr3, rho-rho[k])
                plt.plot(zr[k], 0, 'o')
                plt.plot(zstar(rho), rho-rho[k])
                plt.axis([z0[k]-300, z0[k]+300, -.5, .5])

        zrr = interp1d(r0, z0, bounds_error=False,
                       kind='cubic', fill_value='extrapolate')

        zrho = zrr(rho)
        dzstar = (zrho-z0)
        drho = rho-r0
        dz2 = dzstar/(1.+r0*g*K*dzstar/drho)
        eape0 = 0.5*dz2*drho

        plt.figure(1)
        plt.clf()
        plt.subplot(221)
        plt.plot(r0, -z0, lw=4)
        plt.plot(rho, -z0)
        plt.plot(r0star, -z0)
        plt.plot(rhostar, -z0)
        plt.title('%i - %i' % (i0, i1))

        plt.subplot(223)
        plt.plot(zr3-z0, -z0, 'o-')
    #    plt.plot(dz2, -z0)
    #    plt.plot(zr2-z0, -z0, '+-')
        plt.plot(zr-z0, -z0, '+')
        plt.xlabel('dz')

        plt.subplot(224)
        # plt.plot(rho-r0, -z0)
        # plt.xlabel('drho')
        plt.plot(cs[idx], -z0)
        plt.plot(cs1[idx], -z0)
        plt.xlabel('cs')

        plt.subplot(222)
        plt.plot(eape3, -z0)
    #    plt.plot(eape0, -z0)
        plt.plot(eape2, -z0, '+')
        plt.xlabel('EAPE')

        print('EAPE > 0 is ', np.all(eape[wherenotnan(eape)] >= 0))

        plt.figure(2)
        plt.clf()
        rr = zrefprofiles[b'RHO'][idxpos, :]
        zr4, eape4 = compute_eape(zref, rhoref, cs, rr[:, :])
        plt.plot(zr4.T-zref[:, np.newaxis], -zref, 'k', alpha=0.2)
        plt.plot(np.nanmean(zr4.T-zref[:, np.newaxis], axis=1), -zref, 'r')
