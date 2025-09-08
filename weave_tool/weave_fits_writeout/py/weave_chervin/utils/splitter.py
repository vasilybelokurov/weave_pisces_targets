import h5py
import os
import astropy.table as atpy
import numpy as np
import astropy.io.fits as pyfits


def get_path(dir, fpref, hpx):
    u1 = hpx % 10
    u2 = hpx % 100
    dirpath = f'{dir}/{u1}/{u2}'
    os.makedirs(dirpath, exist_ok=True)
    path = f'{dirpath}/{fpref}_{hpx}.fits'
    return path


def doit(h5file, dir, fpref):
    fp = h5py.File(h5file, 'r')
    nside = 64
    nside0 = 128
    hpx0 = fp['hpx128'][()]
    assert (nside0 >= nside)
    assert (nside0 == (nside0 // nside) * nside)
    hpx = hpx0 // (nside0 // nside)**2
    uhpx = np.unique(hpx)
    assert (np.all(np.diff(hpx) >= 0))
    cols = fp.keys()
    datasets = {}
    for k in cols:
        datasets[k] = fp[k]

    for u in uhpx:
        print('Doing', u)
        # xind = np.nonzero(hpx == u)
        i1 = np.searchsorted(hpx, u, 'left')
        i2 = np.searchsorted(hpx, u, 'right')
        xind = slice(i1, i2)
        path = get_path(dir, fpref, u)
        D = dict([(k, fp[k][xind]) for k in cols])
        T = atpy.Table(D)
        T = pyfits.BinTableHDU(T)
        prim = pyfits.PrimaryHDU()
        prim.header['NSIDE'] = nside
        prim.header['HPX'] = u
        pyfits.HDUList([prim, T]).writeto(path, overwrite=True)
