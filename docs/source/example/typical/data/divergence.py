import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname):
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    return x, y

if __name__ == "__main__":
    argv = sys.argv
    assert 3 == len(argv)
    ifname = argv[1]
    ofname = argv[2]
    x, y = load(ifname)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, color="#FF0000")
    kwrds = {
            "xlabel": "Time",
            "ylabel": "Maximum divergence",
            "ylim": [1.e-16, 1.e-13],
            "yticks": [1.e-16, 1.e-15, 1.e-14, 1.e-13],
            "yscale": "log",
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()

