import astropy.table as atpy
import astropy.coordinates as acoo
import astropy.units as auni
import healpy
import numpy as np

nside0 = 1024  # 0.05 side
nside_weave = 16  # 3.6 side
weave_rad = 1  # deg
pad = 0.2  # fractional padding wrt to radius


class FootPrint:
    def __init__(self):
        self.hpx = np.array([], dtype=int)
        self.nside = nside0

    def add_field(self, ra, dec, rad):
        vec = healpy.ang2vec(ra, dec, lonlat=True)
        pix = healpy.query_disc(nside_weave,
                                vec,
                                np.deg2rad(rad + weave_rad),
                                inclusive=True,
                                nest=True)
        assert self.nside > nside_weave
        nrat = self.nside // nside_weave

        pix = pix[:, None] * nrat**2 + np.arange(nrat**2, dtype=int)[None, :]
        pix = pix.flatten()
        xra, xdec = healpy.pix2ang(self.nside, pix, lonlat=True, nest=True)
        dist = acoo.SkyCoord(ra=xra, dec=xdec, unit=['deg', 'deg']).separation(
            acoo.SkyCoord(ra=ra, dec=dec, unit=['deg', 'deg']))
        xind = dist < (rad + pad) * auni.deg
        self.hpx = np.unique(np.concatenate((self.hpx, pix[xind])))

    def overlaps(self, nside, hpx):
        """
        Check if the given hpx/nside overlaps with the footprint
        """
        assert (nside < self.nside)
        # it only works if I'm querying bigger pixel
        rat = (self.nside // nside)**2
        i1 = np.searchsorted(self.hpx, hpx * rat, 'left') - 1
        i2 = np.searchsorted(self.hpx, hpx * rat + rat - 1, 'right') - 1
        if i1 != i2:
            return True

    def contains_pos(self, ra, dec):
        """
        Check if the given ra,dec is contained in the footprint
        """
        assert (np.isfinite(ra + dec).all())
        hpx = healpy.ang2pix(self.nside, ra, dec, lonlat=True, nest=True)
        return np.isin(hpx, self.hpx)


def read_fields(fits_file):
    tab = atpy.Table().read(fits_file)
    foot = FootPrint()
    for i in range(len(tab)):
        foot.add_field(tab['RA'][i], tab['DEC'][i], tab['RADIUS'][i])
    return foot
