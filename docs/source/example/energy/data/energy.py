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
    for fname in fnames:
        data = np.loadtxt(fname)
        t  = data[:, 0]
        ex = data[:, 1]
        ey = data[:, 2]
        ez = data[:, 3]
        ref = ex[0] + ey[0] + ez[0]
        ls.append(fname.strip().split("energy-")[1].split(".dat")[0])
        xs.append(t[1:])
        ys.append(ex[1:] + ey[1:] + ez[1:] - ref)
    return ls, xs, ys

if __name__ == "__main__":
    argv = sys.argv
    assert(4 == len(argv))
    idname  = argv[1]
    ofname1 = argv[2]
    ofname2 = argv[3]
    ls, xs, ys = load(idname)
    #
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    for l, x, y in zip(ls, xs, ys):
        ax.plot(x, -y, label=l)
    kwrds = {
            "title": "",
            "xlabel": "time",
            "ylabel": "$-\Delta E$",
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
    zs = list()
    for l, y in zip(ls, ys):
        xs.append(float(l))
        zs.append(-y[-1])
    xs = np.array(xs)
    zs = np.array(zs)
    ax.plot(xs, zs, marker="o", label="result", color="#FF0000")
    ax.plot(xs, 2. * zs[0] / xs[0]**3. * xs**3., marker="+", label="3rd order", color="#0000FF")
    kwrds = {
            "title": "",
            "xlabel": "$\Delta t$",
            "ylabel": "$-\Delta E$",
            "xscale": "log",
            "yscale": "log",
    }
    ax.legend()
    ax.set(**kwrds)
    pyplot.savefig(ofname2)
    pyplot.close()

