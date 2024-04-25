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
    ri = 1.
    ui = 1.
    re = np.load(f"{config}/Re.npy")
    lengths = np.load(f"{config}/lengths.npy")
    ro = ri + lengths[0]
    factor = 2. / re * ri**2. * ro**2. / (ro**2. - ri**2.) * (ui / ri)**2. * lengths[1] * lengths[2]
    return values / factor

if __name__ == "__main__":
    argv = sys.argv
    assert 5 == len(argv)
    config = argv[1]
    ifnames = [argv[2], argv[3]]
    ofname = argv[4]
    x0, y0 = load(ifnames[0])
    x1, y1 = load(ifnames[1])
    y0 = normalise(config, y0)
    y1 = normalise(config, y1)
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x0, y0, color="#FF0000", label="dissipation")
    ax.plot(x1, y1, color="#0000FF", label="injection")
    ax.plot(x0, np.full(x0.shape, 1.1375), color="#000000", linestyle="--", label="Ostilla et al., JFM, 2013")
    kwrds = {
            "xlabel": "Time",
            "ylabel": "Nu",
            "xlim": [0., 200.],
    }
    ax.set(**kwrds)
    ax.legend()
    pyplot.savefig(ofname)
    pyplot.close()

