import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PPoly, lagrange


class PLagrange(PPoly):
    def __init__(self, x, y, order=2, extralims=None):
        nx = len(x)+1
        xe = np.zeros((nx+1,))
        xe[1:nx] = x.copy()
        if extralims is None:
            xe[0] = 2*x[0]-x[1]
            xe[-1] = 2*x[-1]-x[-2]
        else:
            xe[0] = extralims[0]
            xe[-1] = extralims[1]
        coef = np.zeros((order, nx, order))
        shift = order//2
        yi = np.zeros((order,))
        for k in range(nx):
            i0 = k-shift
            i1 = i0+order
            if i0 < 0:
                i0 = 0
            if i1 > nx-1:
                i1 = nx-1
            localorder = i1-i0
            for i in range(i0, i1):
                yi *= 0
                yi[i-i0] = y[i]
                p = lagrange(x[i0:i1]-xe[k], yi)
                # print(k, i, xi, yi, p.coef)
                # print(p.coef, localorder)
                # print(len(p.coef), k, i-i0, np.shape(coef))
                coef[order-localorder:, k, i-i0] = p.coef

        super(PLagrange, self).__init__(coef, xe, extrapolate=False)

    def __call__(self, x):
        out = super(PLagrange, self).__call__(x)
        return np.sum(out, axis=-1)


def f(x):
    return np.cos(x/1.5)*np.exp(-x/5)


if __name__ == '__main__':

    xi = np.linspace(0, 10, 7)
    yi = (xi-5)**2/25.  #
    yi = f(xi)

    g = PLagrange(xi, yi, extralims=[-1, 12], order=6)

    x = np.linspace(-1, 12, 51)

    plt.clf()
    plt.plot(x, f(x))
    plt.plot(x, g(x), '.')
    plt.plot(xi, yi, '+', markersize=10)
