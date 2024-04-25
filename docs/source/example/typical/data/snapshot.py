import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
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
    z = np.linspace(0.5 * dz, lz - 0.5 * dz, glksize)
    ux = np.average(ux, axis=1)
    uy = np.average(uy, axis=1)
    uz = np.average(uz, axis=1)
    ux = 0.5 * ux[:, 1:] + 0.5 * ux[:, :-1]
    uz = uz[:, 1:-1]
    return x, z, ux, uy, uz

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    idname = argv[1]
    ofname = argv[2]
    x, z, ux, uy, uz = load(idname)
    fig = pyplot.figure(figsize=(8, 4))
    ax121 = fig.add_subplot(121)
    ax122 = fig.add_subplot(122)
    ax121.contourf(x, z, uy, vmin=0., vmax=1., levels=11, cmap="bwr")
    ax122.contourf(x[1:-1], z, (ux**2. + uz**2.)**0.5, levels=11, cmap="bwr")
    ax122.quiver(x[1:-1], z, ux, uz)
    kwrds = {
            "title": "",
            "aspect": "equal",
            "xlim": [1., 2.],
            "ylim": [0., 2.],
            "xlabel": "r",
            "ylabel": "z",
            "xticks": [1., 1.5, 2.],
            "yticks": [0., 1., 2.],
    }
    ax121.set(**kwrds)
    ax122.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()

