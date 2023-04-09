import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def load(fname, xmin):
    data = np.loadtxt(fname)
    x  = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    y0 = y0[x > xmin]
    y1 = y1[x > xmin]
    x  =  x[x > xmin]
    return x, y0, y1

if __name__ == "__main__":
    argv = sys.argv
    assert(3 == len(argv))
    ifname = argv[1]
    ofname = argv[2]
    xmin = 100.
    x, y0, y1 = load(ifname, xmin)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y0, color="#FF0000", label="inner cylinder")
    ax.plot(x, y1, color="#0000FF", label="outer cylinder")
    # Ostilla et al., JFM, 2013
    y2 = np.full(x.shape, 1.1375)
    ax.plot(x, y2, color="#000000", linestyle="--", label="Ostilla et al., JFM, 2013")
    kwrds = {
            "title": "",
            "xlim": [xmin, 200.],
            "xlabel": "time",
            "ylabel": "Normalised torque",
            "xticks": [100., 150., 200.],
    }
    ax.set(**kwrds)
    ax.legend()
    pyplot.savefig(ofname)
    pyplot.close()

