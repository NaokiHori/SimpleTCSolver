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

def normalise(config, values):
    return values / factor

if __name__ == "__main__":
    argv = sys.argv
    assert 4 == len(argv)
    ifnames = [argv[1], argv[2]]
    ofname = argv[3]
    x0, y0 = load(ifnames[0])
    x1, y1 = load(ifnames[1])
    y = np.abs(y1 - y0)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x0, y, color="#FF0000")
    kwrds = {
            "title": "",
            "xlabel": "Time",
            "ylabel": "Injection - Dissipation",
            "yscale": "log",
    }
    ax.set(**kwrds)
    pyplot.savefig(ofname)
    pyplot.close()

