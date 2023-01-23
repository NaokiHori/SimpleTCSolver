# script used to draw bottom picture in README

import sys
import numpy as np
from matplotlib import pyplot as plt


def flowfield(root, ax):
    xc = np.load(f"{root}/xc.npy")
    zc = np.load(f"{root}/zc.npy")
    xf = np.load(f"{root}/xf.npy")
    zf = np.load(f"{root}/zf.npy")
    uy = np.load(f"{root}/uy.npy")
    # average in theta
    uy = np.average(uy, axis=1)
    ax.contourf(xc, zc, uy, vmin=0., vmax=1., levels=11, cmap="bwr")
    keywords = {
            "title": "Azimuthal velocity",
            "aspect": "equal",
            "xlim": [xf[0], xf[-1]],
            "ylim": [zf[0], zf[-1]],
            "xlabel": "",
            "ylabel": "",
            "xticks": [],
            "yticks": [],
    }
    ax.set(**keywords)

def torque(root, ax):
    data = np.loadtxt(f"{root}/nusselt.dat")
    t_ = data[:, 0]
    v0 = data[:, 1]
    v1 = data[:, 2]
    t  = t_[t_ > 100.]
    v0 = v0[t_ > 100.]
    v1 = v1[t_ > 100.]
    vr = np.full(t.shape, 1.1375)
    ax.set_title("Normalised torque")
    ax.plot(t, v0, linestyle="solid",  label="inner")
    ax.plot(t, v1, linestyle="solid",  label="outer")
    ax.plot(t, vr, linestyle="dashed", label="Ostilla et al., JFM, 2013")
    ax.legend()

if __name__ == "__main__":
    argv = sys.argv
    assert(len(argv) == 3)
    root = argv[1]
    filename = argv[2]
    fig = plt.figure(figsize=(8, 6))
    ax121 = fig.add_subplot(121)
    ax122 = fig.add_subplot(122)
    flowfield(f"{root}/flowfield", ax121)
    torque(root, ax122)
    plt.savefig(filename, dpi=300)
    plt.close()

