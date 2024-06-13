import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(dname):
    fnames = [f"{dname}/{fname}" for fname in os.listdir(dname) if fname.startswith("energy") and fname.endswith(".dat")]
    fnames = sorted(fnames)
    ls = list()
    xs = list()
    ys = list()
    zs = list()
    for fname in fnames:
        data = np.loadtxt(fname)
        t  = data[:, 0]
        ex = data[:, 1]
        ey = data[:, 2]
        ez = data[:, 3]
        sc = data[:, 4]
        et = ex + ey + ez
        e0 = et[0]
        s0 = sc[0]
        et /= e0
        sc /= s0
        ls.append(fname.strip().split("energy-")[1].split(".dat")[0])
        xs.append(t[1:])
        ys.append(1. - et[1:])
        zs.append(1. - sc[1:])
    return ls, xs, ys, zs

if __name__ == "__main__":
    argv = sys.argv
    assert(4 == len(argv))
    idname  = argv[1]
    ofname1 = argv[2]
    ofname2 = argv[3]
    ls, xs, ys, zs = load(idname)
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    colors = pyplot.rcParams["axes.prop_cycle"].by_key()["color"]
    for l, x, y, z, color in zip(ls, xs, ys, zs, colors):
        ax.plot(x, y, linestyle="solid",  color=color, label=l)
        ax.plot(x, z, linestyle="dashed", color=color)
    kwrds = {
            "title": "",
            "xlabel": "time",
            "ylabel": "$-\Delta E / E_0$",
            "yscale": "log",
    }
    ax.legend()
    ax.set(**kwrds)
    pyplot.savefig(ofname1)
    pyplot.close()
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    xs = list()
    y1s = list()
    z1s = list()
    for l, y, z in zip(ls, ys, zs):
        xs.append(float(l))
        y1s.append(y[-1])
        z1s.append(z[-1])
    xs = np.array(xs)
    ax.plot(xs, y1s, marker="o", label="velocity", color="#FF0000")
    ax.plot(xs, z1s, marker="o", label="scalar",   color="#0000FF")
    ax.plot(xs, 2. * y1s[0] / xs[0]**3. * xs**3., label="3rd order", color="#000000", linestyle="dashed")
    kwrds = {
            "title": "",
            "xlabel": "$\Delta t$",
            "ylabel": "$-\Delta E / E_0$",
            "xscale": "log",
            "yscale": "log",
    }
    ax.legend()
    ax.set(**kwrds)
    pyplot.savefig(ofname2)
    pyplot.close()

