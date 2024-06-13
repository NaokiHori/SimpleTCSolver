import sys
import numpy as np
import matplotlib
# matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname):
    x = np.load(f"{dname}/xc.npy")
    lengths = np.load(f"{dname}/lengths.npy")
    glsizes = np.load(f"{dname}/glsizes.npy")
    lz = lengths[2]
    glksize = glsizes[2]
    dz = lz / glksize
    ux = np.load(f"{dname}/ux.npy")
    uy = np.load(f"{dname}/uy.npy")
    uz = np.load(f"{dname}/uz.npy")
    t = np.load(f"{dname}/t.npy")
    z = np.linspace(0.5 * dz, lz - 0.5 * dz, glksize)
    ux = np.average(ux, axis=1)
    uy = np.average(uy, axis=1)
    uz = np.average(uz, axis=1)
    t = np.average(t, axis=1)
    ux = 0.5 * ux[:, 1:] + 0.5 * ux[:, :-1]
    uz = uz[:, 1:-1]
    # enforce plumes to be centre
    argmin = np.argmin(np.average(uy, axis=1))
    ux = np.roll(ux, shift=-argmin, axis=0)
    uy = np.roll(uy, shift=-argmin, axis=0)
    uz = np.roll(uz, shift=-argmin, axis=0)
    t = np.roll(t, shift=-argmin, axis=0)
    return x, z, ux, uy, uz, t

if __name__ == "__main__":
    argv = sys.argv
    assert 3 == len(argv)
    idname = argv[1]
    ofname = argv[2]
    x, z, ux, uy, uz, t = load(idname)
    fig = pyplot.figure()
    ax131 = fig.add_subplot(131)
    ax132 = fig.add_subplot(132)
    ax133 = fig.add_subplot(133)
    ax131.contourf(x, z, uy, vmin=0., vmax=1., levels=11, cmap="bwr")
    ax132.contourf(x[1:-1], z, (ux**2. + uz**2.)**0.5, levels=11, cmap="bwr")
    ax132.quiver(x[1:-1], z, ux, uz)
    ax133.contourf(x, z, t, vmin=-0.5, vmax=+0.5, levels=11, cmap="viridis")
    kwrds = {
            "aspect": "equal",
            "xlim": [1., 2.],
            "ylim": [0., 2.],
            "xlabel": "r",
            "ylabel": "z",
            "xticks": [1., 1.5, 2.],
            "yticks": [0., 1., 2.],
    }
    ax131.set(**kwrds)
    ax132.set(**kwrds)
    ax133.set(**kwrds)
    ax131.set_title("uy")
    ax132.set_title("roll")
    ax133.set_title("scalar")
    pyplot.savefig(ofname)
    pyplot.close()

