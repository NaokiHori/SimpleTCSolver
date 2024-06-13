import os
import sys
import numpy as np


rng = np.random.default_rng()


def get_is_curved():
    truish_list = ["True", "true"]
    is_curved = os.environ.get("is_curved")
    if None == is_curved:
        print("is_curved is not given")
        sys.exit(1)
    if is_curved in truish_list:
        return True
    else:
        return False


def get_lengths():
    lx = os.environ.get("lx")
    ly = os.environ.get("ly")
    lz = os.environ.get("lz")
    if lz:
        lengths = [lx, ly, lz]
    else:
        lengths = [lx, ly]
    return list(map(float, lengths))


def get_glsizes():
    glisize = os.environ.get("glisize")
    gljsize = os.environ.get("gljsize")
    glksize = os.environ.get("glksize")
    if glksize:
        glsizes = [glisize, gljsize, glksize]
    else:
        glsizes = [glisize, gljsize]
    return list(map(int, glsizes))


def init_time(dest):
    # iterator and time
    step = np.array(0, dtype=np.uint64)
    time = np.array(0, dtype=np.float64)
    np.save(f"{dest}/step.npy", step)
    np.save(f"{dest}/time.npy", time)
    return


def init_domain(is_curved, lengths, glsizes, dest):
    # generate equidistant sequence
    # NOTE: cell face has +1 elements
    xf = np.arange(0, glsizes[0] + 1, 1)
    # stretched grid, clipped Chebyshev just as an example
    is_uniform = False
    if not is_uniform:
        # number of grid points to be clipped at the edges
        nclip = 3
        # gather close to the boundaries
        xf = np.cos(np.pi * (xf + 1. * nclip) / (glsizes[0] + 2. * nclip))
        # make the descending order ascending
        xf *= -1.
    # normalse to enforce [0 : lx]
    xf = lengths[0] * (xf - np.min(xf)) / (np.max(xf) - np.min(xf))
    # cell centers are located at the center
    #   of the two neighbouring cell faces,
    #   which are appended by the boundaries
    xc = 0.
    xc = np.append(xc, 0.5 * xf[:-1] + 0.5 * xf[1:])
    xc = np.append(xc, lengths[0])
    # offset
    # inner cylinder radius for TC domains, otherwise (planar channel) just set 0
    xmin = 1. if is_curved else 0.
    xf += xmin
    xc += xmin
    np.save(f"{dest}/is_curved.npy", np.array(is_curved, dtype=np.bool_))
    np.save(f"{dest}/xf.npy", np.array(xf, dtype=np.float64))
    np.save(f"{dest}/xc.npy", np.array(xc, dtype=np.float64))
    np.save(f"{dest}/glsizes.npy", np.array(glsizes, dtype=np.uint64))
    np.save(f"{dest}/lengths.npy", np.array(lengths, dtype=np.float64))
    return xf, xc


def init_fluid(is_3d, lengths, glsizes, xf, xc, dest):
    if is_3d:
        ux = rng.random((glsizes[2], glsizes[1], glsizes[0] + 1), dtype=np.float64)
        uy = rng.random((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
        uz = rng.random((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
        p  = rng.random((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
        t  = rng.random((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
        ux -= 0.5
        uz -= 0.5
        p  -= 0.5
        t  -= 0.5
        np.save(f"{dest}/ux.npy", ux)
        np.save(f"{dest}/uy.npy", uy)
        np.save(f"{dest}/uz.npy", uz)
        np.save(f"{dest}/p.npy", p)
        np.save(f"{dest}/t.npy", t)
    else:
        ux = rng.random((glsizes[1], glsizes[0] + 1), dtype=np.float64)
        uy = rng.random((glsizes[1], glsizes[0] + 2), dtype=np.float64)
        p  = rng.random((glsizes[1], glsizes[0] + 2), dtype=np.float64)
        t  = rng.random((glsizes[1], glsizes[0] + 2), dtype=np.float64)
        ux -= 0.5
        p  -= 0.5
        t  -= 0.5
        np.save(f"{dest}/ux.npy", ux)
        np.save(f"{dest}/uy.npy", uy)
        np.save(f"{dest}/p.npy", p)
        np.save(f"{dest}/t.npy", t)


def main():
    is_curved = get_is_curved()
    lengths = get_lengths()
    glsizes = get_glsizes()
    assert len(lengths) == len(glsizes)
    dest = sys.argv[1]
    is_3d = 3 == len(lengths)
    if is_3d:
        print("A 3D field is initialised")
    else:
        print("A 2D field is initialised")
    # init and save
    init_time(dest)
    xf, xc = init_domain(is_curved, lengths, glsizes, dest)
    init_fluid(is_3d, lengths, glsizes, xf, xc, dest)


main()
