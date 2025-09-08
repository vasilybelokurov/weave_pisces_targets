import os
import glob
import math
import numpy as np
import healpy
import scipy
import yaml
import scipy.stats
import astropy.table as atpy
import astropy.io.fits as pyfits
import astropy.coordinates as acoo
import astropy.units as auni
import weave_galr.utils.footprint
import fitsio
import dustmaps.sfd
from . import gaia_extinction

acoo.galactocentric_frame_defaults.set('v4.0')
# This is to fix the definition of GC frame

SFDQuery = dustmaps.sfd.SFDQuery()


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


weave_fov = 3.14
weave_rad = np.sqrt(weave_fov / np.pi)


def getParams(fname=None):
    if fname is None:
        confdir = os.path.dirname(__file__) + '/../conf/'
        fname = confdir + '/selection.yaml'
    with open(fname, 'r') as fp:
        SelConfig = yaml.load(fp, Loader=yaml.FullLoader)
    SelParam = SelConfig['SelectionParameters']
    SelParamSV = SelConfig['SelectionParametersSV']

    return SelParam, SelParamSV


class Mask:
    def __init__(self):
        self.ras = []
        self.decs = []
        self.rads = []

    def add(self, ra, dec, rad):
        self.ras.append(ra)
        self.decs.append(dec)
        self.rads.append(rad)

    def inside(self, ra, dec):
        Cmask = acoo.SkyCoord(ra=np.array(self.ras) * auni.deg,
                              dec=np.array(self.decs) * auni.deg)
        Cdat = acoo.SkyCoord(ra=np.array(ra) * auni.deg,
                             dec=np.array(dec) * auni.deg)
        sep = Cdat[:, None].separation(Cmask[None, :])
        ret_inside = np.any(
            sep.to_value(auni.deg) < np.array(self.rads)[None, :], axis=1)
        return ret_inside


def load_mask():
    fname = os.path.dirname(__file__) + '/../data/SharedSelectionMask.csv'
    T = atpy.Table().read(fname)
    SkyMask = Mask()
    for i in range(len(T)):
        SkyMask.add(T['ra'][i], T['dec'][i], T['radius'][i])
    return SkyMask


SkyMask = load_mask()


def get_selconfig_dict(SelParam, SelParamSV):
    res = {}
    paramlist = [('SelectionParam', SelParam),
                 ('selectionParamSV', SelParamSV)]
    for name, param in paramlist:
        for k, v in param.items():
            if not isinstance(v, dict):
                res[name + '_' + k] = v
            else:
                for k1, v1 in v.items():
                    res[name + '_' + k + '_' + k1] = v1
    return res


def fill_selconfig_header(hdr, SelParam, SelParamSV):
    for k, v in get_selconfig_dict(SelParam, SelParamSV).items():
        hdr['HIERARCH ' + k] = v


def good_object(tab):
    star_galaxy = ((tab['rmeanpsfmag'] - tab['rmeankronmag']) < 0.05)
    reasonable_color = ((tab['gmeanpsfmag'] - tab['rmeanpsfmag']) > -1) & (
        (tab['gmeanpsfmag'] - tab['rmeanpsfmag']) < 4)
    aen_sel = (tab['astrometric_excess_noise'] <=
               1) | (~np.isfinite(tab['astrometric_excess_noise']))
    return star_galaxy & reasonable_color & aen_sel


def get_rand_subset(mask, n, rstate=None):
    # get a random subset of n elements out masked set
    xids = np.nonzero(mask)[0]
    retmask = np.zeros_like(mask)
    xids = rstate.permutation(xids)
    retmask[xids[:n]] = True
    return retmask


def get_pm_ref(ra, dec, dist):
    """ get the reference proper motion of a non-rotating population
    at distance dist kpc and ra,dec
    """
    coo = acoo.SkyCoord(ra=ra * auni.deg,
                        dec=dec * auni.deg,
                        distance=dist * auni.kpc)
    coo1 = coo.transform_to(acoo.Galactocentric())
    kms = auni.km / auni.s
    coo2 = acoo.Galactocentric(x=coo1.x,
                               y=coo1.y,
                               z=coo1.z,
                               v_x=0 * kms,
                               v_y=0 * kms,
                               v_z=0 * kms)
    coo3 = coo2.transform_to(acoo.ICRS())
    return (coo3.pm_ra_cosdec.to_value(auni.mas / auni.year),
            coo3.pm_dec.to_value(auni.mas / auni.year))


def get_Gabs_array():
    X = np.array(
        [[
            0.35295079, 0.35837824, 0.36494692, 0.37081186, 0.37644827,
            0.3823219, 0.38850216, 0.39421875, 0.39991164, 0.40581852,
            0.41160093, 0.41742374, 0.4232192, 0.42884001, 0.43507031,
            0.44088076, 0.44665666, 0.45243954, 0.45831627, 0.46420393,
            0.47008122, 0.47591683, 0.48164511, 0.48784672, 0.49335117,
            0.49970303, 0.50554446, 0.51128497, 0.51712885, 0.52304527,
            0.52865908, 0.53482021, 0.54079635, 0.54628985, 0.55209782,
            0.55761943, 0.56384889, 0.56964285, 0.57588722, 0.58166575,
            0.58731464, 0.59289463, 0.59879211, 0.60506157, 0.61067938,
            0.61645146, 0.62320408, 0.62864145, 0.63398677, 0.63980294,
            0.64581897, 0.65172987, 0.65717409, 0.66363244, 0.66959398,
            0.67502135, 0.68124196, 0.68718534, 0.6925811, 0.69830519,
            0.70437427, 0.71067416, 0.71585453, 0.72163193, 0.72852886,
            0.73391497, 0.73916881, 0.74546199, 0.75117916, 0.75731704,
            0.76298121, 0.76885006, 0.77505648, 0.7807736, 0.78595958,
            0.79255705, 0.79850522, 0.80458335, 0.8100115, 0.81576132,
            0.82169893, 0.82765842, 0.83279645, 0.83948452, 0.84540511,
            0.85063794, 0.85694483, 0.86270215, 0.86848022, 0.87406537,
            0.88015075, 0.88614147, 0.89197837, 0.89805024, 0.90373473,
            0.90951077, 0.91571302, 0.92137301, 0.92727303, 0.93306249,
            0.93903272, 0.94506703, 0.9503768, 0.95633258, 0.9624521,
            0.96783988, 0.97473933, 0.97972362, 0.98497591, 0.99196318,
            0.9979141, 1.00286428, 1.0090607, 1.0142569, 1.02062108,
            1.02698995, 1.03249851, 1.03851082, 1.04344773, 1.04979027,
            1.05597308, 1.06137017, 1.06768043, 1.07450218, 1.07977691,
            1.08540328, 1.08998587, 1.09789533
        ],
         [
             4.38051582, 4.20826533, 4.23208934, 4.31437742, 4.229072,
             4.24233475, 4.2521481, 4.34353879, 4.36275193, 4.38980165,
             4.42672319, 4.52376109, 4.62198845, 4.69206753, 4.7432576,
             4.7907839, 4.88973865, 5.0861731, 5.20021283, 5.23329925,
             5.30240503, 5.38894607, 5.44811237, 5.53574141, 5.6285056,
             5.65494247, 5.70257369, 5.76503454, 5.79145142, 5.88658943,
             5.94132958, 5.99191386, 6.06435782, 6.08790393, 6.12264819,
             6.20486633, 6.23985609, 6.26786472, 6.35044378, 6.38346919,
             6.46962684, 6.44758747, 6.54489066, 6.56497229, 6.59635288,
             6.68601488, 6.73969697, 6.7823446, 6.82257707, 6.8517389,
             6.89395448, 6.96881779, 6.98724052, 7.02961238, 7.0707874,
             7.14949638, 7.1507643, 7.19416196, 7.23652777, 7.29061699,
             7.36890584, 7.3737611, 7.4776653, 7.55710189, 7.58648036,
             7.57940612, 7.59528912, 7.75936793, 7.79808996, 7.82290071,
             7.9238059, 7.84242117, 8.00528331, 8.00993803, 8.17008743,
             8.20886189, 8.18099346, 8.30319152, 8.36350428, 8.3836347,
             8.58274501, 8.50110759, 8.63750526, 8.65857414, 8.75836527,
             8.79240947, 8.91287532, 9.0026598, 9.1474261, 9.07368256,
             9.12464961, 9.2797731, 9.26842964, 9.36374644, 9.5022904,
             9.54556295, 9.57518125, 9.61475196, 9.63988651, 9.87524017,
             9.82569435, 9.98292033, 10.16558703, 9.9444265, 10.28861896,
             10.37511487, 10.36858534, 10.39457935, 10.36511662, 10.56770023,
             10.66896779, 10.68227931, 10.90067812, 10.93779758, 10.95675811,
             10.95964365, 11.16238105, 11.16651184, 11.25826073, 11.44707216,
             11.46920113, 11.64274656, 11.4047966, 11.73148317, 11.69876544,
             11.80029814, 11.87049207, 11.82816339
         ],
         [
             0.16347342, 0.21071131, 0.18063469, 0.25291923, 0.19937065,
             0.16890873, 0.18127217, 0.21980206, 0.24325483, 0.25007444,
             0.26041935, 0.26622234, 0.28826367, 0.29905599, 0.33191121,
             0.35581726, 0.34694151, 0.24130378, 0.16307402, 0.18219049,
             0.15037215, 0.17631211, 0.17382735, 0.17363072, 0.19170819,
             0.14923396, 0.15602657, 0.14371296, 0.17628815, 0.14820106,
             0.14773153, 0.12941767, 0.14842124, 0.16256556, 0.13878892,
             0.13192525, 0.13487455, 0.13846747, 0.13299956, 0.13723199,
             0.15911344, 0.14863377, 0.14749634, 0.13591851, 0.11785963,
             0.14569317, 0.14584796, 0.13984712, 0.12443507, 0.14289908,
             0.12226686, 0.14161121, 0.17161486, 0.14230318, 0.16551339,
             0.14658724, 0.12418722, 0.13578857, 0.11869855, 0.15092282,
             0.16734492, 0.1685851, 0.2845134, 0.21639416, 0.31590682,
             0.22259608, 0.15235285, 0.27601321, 0.3807229, 0.27021635,
             0.26741, 0.28111641, 0.36115087, 0.28624189, 0.34210404,
             0.31296406, 0.38290126, 0.42823856, 0.30402349, 0.38406077,
             0.44066601, 0.3800203, 0.43489793, 0.38464327, 0.44126347,
             0.3734938, 0.45450335, 0.49403045, 0.51159537, 0.5132253,
             0.4859085, 0.6025285, 0.45483643, 0.5190396, 0.56110865,
             0.42978032, 0.53651409, 0.48228044, 0.48163477, 0.63606138,
             0.54343266, 0.46294575, 0.51791147, 0.45836082, 0.37308321,
             0.53084576, 0.44635038, 0.56621238, 0.41478257, 0.35788351,
             0.40068817, 0.38262741, 0.44392008, 0.39567206, 0.4905301,
             0.38990816, 0.39215462, 0.32650903, 0.39816253, 0.33870536,
             0.35547549, 0.38953679, 0.17670484, 0.33636289, 0.30609868,
             0.40181109, 0.08762758, 0.14434961
         ]])
    return X


def gaiapoe_sel(**kw):
    SelParam = kw['SelParam']
    SelParamSV = kw['SelParamSV']
    parallax = kw['parallax']
    POE = parallax / kw['parallax_error']
    zpt = 17. / 1000.
    dist = 1 / (parallax + zpt)
    AG = gaia_extinction.get(kw['phot_g_mean_mag'],
                             kw['phot_bp_mean_mag'],
                             kw['phot_rp_mean_mag'],
                             kw['ebv'],
                             maxnit=1)[0]
    # Absolute magnitude needs extinction in the G-band "AG"
    dist[dist <= 0] = np.nan
    gmag = kw['phot_g_mean_mag']
    MG = gmag - 5 * np.log10(1000 * dist) + 5 - AG
    Vtan = dist * np.sqrt(kw['pmra']**2 + kw['pmdec']**2) * 4.74047
    minmag = SelParam['POE']['min_gaia_mag']
    maxmag = SelParam['POE']['max_gaia_mag']
    mask = ((dist < 3.0) & (MG > 3.2) & (POE >= 5.0) & betw(Vtan, 200, 800) &
            (kw['ruwe'] < 1.4)) & betw(gmag, minmag, maxmag)
    Ntarget = kw['Ntarget']
    return get_rand_subset(mask, Ntarget, rstate=kw['rstate'])


def add_Gabs_fit(color, Nbins_fit=128):
    """
    Gets the absolute G-magnitude form G - RP colours.
    NOTE: uses an external file with mean MS locus and dispersion.
    """
    g_rp = color.copy()
    hist, edges = np.histogram(g_rp, bins=Nbins_fit, range=(0.35, 1.1))
    centres = (edges[1::] + edges[:-1]) / 2
    g_rp_mode = centres[hist.argmax()]
    g_rp[g_rp < 0.35] = g_rp_mode
    g_rp[g_rp >= 1.1] = 1.099
    g_rp[~np.isfinite(g_rp)] = g_rp_mode
    # NB: limits should be the same as used to obtain mean_GABS
    g_rp_binned = np.digitize(g_rp, edges) - 1
    g_rp_binned[g_rp_binned < 0] = 0
    # array with Gabs and error for stars in each of the 128 bins
    MSfit = get_Gabs_array()
    MSmean, MSstd = MSfit[[0, 1]], MSfit[2]
    assert (len(MSstd) == Nbins_fit)
    # Returning the fitted Gabs, and with of the MS.
    return np.array(MSmean[1])[g_rp_binned], np.array(MSstd)[g_rp_binned]


def prop_Hg_uncertainty(**kw):
    """
    Propagate uncertainties for Hg, ignores uncertainty in G-magnitude.
    Procedurally generated, apologies for the long expression.
    """
    from numpy import sqrt, log

    pmra = kw['pmra']
    pmdec = kw['pmdec']
    pmra_error = kw['pmra_error']
    pmdec_error = kw['pmdec_error']
    pmra_pmdec_corr = kw['pmra_pmdec_corr']

    return sqrt(
        ((((((5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                   ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmdec)))) *
             (pmdec_error**2)) *
            (5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                  ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmdec))))) +
           (((((5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) * (
               (0.5 * (((pmra**2) + (pmdec**2))**-0.5)) *
               (2 * pmdec)))) * pmra_pmdec_corr) * pmra_error) * pmdec_error) *
            (5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                  ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmra)))))) +
          (((((5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) * (
              (0.5 * (((pmra**2) + (pmdec**2))**-0.5)) *
              (2 * pmra)))) * pmra_pmdec_corr) * pmra_error) * pmdec_error) *
           (5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                 ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmdec)))))) +
         (((5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                 ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmra)))) *
           (pmra_error**2)) *
          (5 * ((1 / (sqrt(((pmra**2) + (pmdec**2))) * log(10))) *
                ((0.5 * (((pmra**2) + (pmdec**2))**-0.5)) * (2 * pmra)))))))


def gaiarpm_sel(**kw):
    '''
    Compute distance using RPM method see Viswanathan et al., in prep and
    Koppelman+2019 (https://arxiv.org/abs/2004.07328)
    '''
    SelParam = kw['SelParam']
    SelParamSV = kw['SelParamSV']

    AG, ABP, ARP = gaia_extinction.get(kw['phot_g_mean_mag'],
                                       kw['phot_bp_mean_mag'],
                                       kw['phot_rp_mean_mag'],
                                       kw['ebv'],
                                       maxnit=1)

    gmag = kw['phot_g_mean_mag']
    G_rp = (gmag - AG) - (kw['phot_rp_mean_mag'] - ARP)
    G = gmag - AG
    Gabs = G + 5 * np.log10(np.maximum(kw['parallax'], 1e-6)) - 10
    PM = np.sqrt(kw['pmra']**2 + kw['pmdec']**2)
    Hg = G + 5 * np.log10(PM) - 10
    Hg_uncertainty = prop_Hg_uncertainty(**kw)

    # Some masks
    POE5 = (kw['parallax'] / kw['parallax_error'] > 5)
    quality = (kw['ruwe'] < 1.4) & (AG < 2) & (Hg > Hg_uncertainty * 10**1.75)

    # three linear fits to make cuts on the RPM diagram for the halo selection
    p1 = [11.61608077, -0.24101326]
    p2 = [8.13255287, 1.6000506]
    p3 = [11.49163297, -0.87464927]

    # the lines
    L1 = p1[1] + p1[0] * G_rp
    L2 = p2[1] + p2[0] * G_rp
    L3 = p3[1] + p3[0] * G_rp

    # g_rp bounds
    g_rp_left = 0.35
    g_rp_right = 1.1
    g_rp_split12 = 0.5285055589970487
    g_rp_split23 = 0.7367195185319125

    # tangential velocity bounds
    VT_upper = 200
    VT_lower = 800

    # removing white dwarfs - see Viswanathan et al., in prep
    WDS1 = (Gabs < L1 + 2) & (G_rp < g_rp_split12)
    WDS2 = (Gabs < L2 + 2) & (G_rp > g_rp_split12) & (G_rp < g_rp_split23)
    WDS3 = (Gabs < L3 + 2) & (G_rp > g_rp_split23)
    WD = ((WDS1 | WDS2 | WDS3) * POE5)

    # selection of halo stars in the RPM diagram
    U1 = (Hg > L1 + 5 * np.log10(VT_upper / 4.74047))
    L1 = (Hg < L1 + 5 * np.log10(VT_lower / 4.74047))
    MSS1 = (G_rp > g_rp_left) & (G_rp < g_rp_split12) & (U1) & (L1)

    U2 = (Hg > L2 + 5 * np.log10(VT_upper / 4.74047))
    L2 = (Hg < L2 + 5 * np.log10(VT_lower / 4.74047))
    MSS2 = (G_rp > g_rp_split12) & (G_rp < g_rp_split23) & (U2) & (L2)

    U3 = (Hg > L3 + 5 * np.log10(VT_upper / 4.74047))
    L3 = (Hg < L3 + 5 * np.log10(VT_lower / 4.74047))
    MSS3 = (G_rp > g_rp_split23) & (G_rp < g_rp_right) & (U3) & (L3)

    MS = ((MSS1 | MSS2 | MSS3) * (quality) * (~WD))

    Grpm, Grpm_std = add_Gabs_fit(G_rp)

    phot_dist = 10**((G - Grpm - 10) / 5)

    # RPM sample uses the same magnitude limites as POE
    minmag = SelParam['POE']['min_gaia_mag']
    maxmag = SelParam['POE']['max_gaia_mag']

    # Return all main-sequence, good astrometry and
    # not WD within a 5kpc volume.
    mask = MS * (phot_dist < 5) * betw(gmag, minmag, maxmag)

    Ntarget = kw['Ntarget']
    return get_rand_subset(mask, Ntarget, rstate=kw['rstate'])


def blue_sel(**kw):
    # function implementing the blue box  selection
    SelParam = kw['SelParam']
    g, r = kw['mag_g'], kw['mag_r'],
    sv = kw.get('sv') or False
    col1, col2 = (SelParam['MSTO']['min_gr_col'],
                  SelParam['MSTO']['max_gr_col'])
    if sv:
        SelParamSV = kw['SelParamSV']
        mag1, mag2 = (SelParam['MSTO']['max_r_mag'],
                      SelParamSV['MSTO']['max_r_mag'])
    else:
        mag1, mag2 = (SelParam['MSTO']['min_r_mag'],
                      SelParam['MSTO']['max_r_mag'])
    ruwe_sel = (~np.isfinite(kw['ruwe'])) | (kw['ruwe'] < 1.4)
    return betw(g - r, col1, col2) & betw(r, mag1, mag2) & ruwe_sel


def filler_sel(**kw):
    # filler selection
    SelParam = kw['SelParam']
    r = kw['mag_r']
    mag1, mag2 = (SelParam['FILLER']['min_r_mag'],
                  SelParam['FILLER']['max_r_mag'])
    return betw(r, mag1, mag2)


def red_sel(**kw):
    # the red box selection
    SelParam = kw['SelParam']
    SelParamSV = kw['SelParamSV']

    # in the case of sv I have extra selection that goes fainter
    minmag, maxmag = (SelParam['GIANTS']['min_r_mag'],
                      SelParam['GIANTS']['max_r_mag'])
    if 'sv' in kw:
        minmag, maxmag = (SelParam['GIANTS']['max_r_mag'],
                          SelParamSV['GIANTS']['max_r_mag'])
    mincol, maxcol = (SelParam['GIANTS']['min_gr_col'],
                      SelParam['GIANTS']['max_gr_col'])

    g, r = kw['mag_g'], kw['mag_r']
    pmra, pmdec = kw['pmra'], kw['pmdec']
    # ra0, dec0 = kw['ra0'], kw['dec0']  # center of the field
    plx, eplx = kw['plx'], kw['eplx']
    pmtot = (pmra**2 + pmdec**2)**.5

    maxvtan = SelParam['GIANTS']['max_vtan']
    mind = SelParam['GIANTS']['min_dist']

    # Fiducial magnitude limit as function of colour
    def MagLim(col):
        return 4 - (col - 0.65) * 6

    # closest distance for star with colour col, magnitude mag (in kpc)
    def DistLim(col, mag, dmin):
        return np.maximum(10**((mag - MagLim(col)) / 5. - 2), dmin)

    def PMErr(mag):
        return 2.5 * 10**(0.27 * (mag - 16) - 1.3)

    curDistLim = DistLim(g - r, r, mind)
    plxlim = 2.5 * eplx + 0.06 + 1. / curDistLim
    pmlim = maxvtan / curDistLim / 4.74 + PMErr(r)

    sel = np.ones(len(g), dtype=bool)
    sel = sel & (kw['ruwe'] < 1.4)
    sel = sel & (plx < plxlim)

    sel = sel & (pmtot < pmlim)
    sel = sel & betw(r, minmag, maxmag)
    sel = sel & betw(g - r, mincol, maxcol)
    return sel


def select_bin(subset0, mag, ntarget, maglim=None, nbins=3, rstate=None):
    '''
    Select stars in a subset, so that they uniformly
    populate the bins in magnitude
    '''

    # limit to stars within
    subset = subset0 & (mag >= maglim[0]) & (mag < maglim[1])
    if subset.sum() < ntarget:
        return subset
    # edges
    bine = np.linspace(maglim[0], maglim[1], nbins + 1)
    hh, loc = np.histogram(mag[subset], bins=bine)
    # positions of stars in the bins
    pos = np.digitize(mag[subset], bine) - 1

    # assigned numbers for a given bins (assuming equal number in each bin)
    # splits = np.zeros(nbins, dtype=int) + ntarget // nbins
    # ensure that the total is exactly ntarget
    # splits[0] = ntarget - (splits[1:].sum())
    splits = rstate.multinomial(ntarget, np.zeros(nbins) + 1. / nbins)
    assert (splits.sum() == ntarget)

    # the ratio of wanted to actual number in each bin
    fracs = splits * 1. / (hh + (hh == 0))
    randoms1 = rstate.uniform(size=subset.sum())

    # this will select all objects in empty bins,
    # and a fraction in overfilled bins
    selected1 = randoms1 < fracs[pos]

    # these are ids within the subset
    xids1 = rstate.permutation(np.nonzero(subset)[0][selected1])
    # these will go first in the list
    xids2 = rstate.permutation(np.nonzero(subset)[0][~selected1])
    # these will go if needed

    # select first ntarget objects
    xids = np.concatenate((xids1, xids2))[:ntarget]

    mask = np.zeros_like(subset, dtype=bool)
    mask[xids] = True
    assert (mask.sum() == ntarget)
    assert ((mask & subset).sum() == ntarget)
    return mask


def main_sel(tab,
             nfibers=600,
             nsv=0,
             rstate=None,
             SelParam=None,
             SelParamSV=None):
    # main function performing selection of the input set of stars
    # I assume that the function takes sa input stars already within
    # the field of view

    select_in_bins = True
    strategy = 'redfrac'
    priority_giants = 10
    ra0, dec0 = [np.median(_) for _ in [tab['ra'], tab['dec']]]
    # field center

    blue_set = blue_sel(
        mag_g=tab['mag_g'],
        mag_r=tab['mag_r'],
        ruwe=tab['ruwe'],
        SelParam=SelParam,
        SelParamSV=SelParamSV,
    )
    # the subset of stars that fits our blue selection

    red_kw = dict(
        mag_g=tab['mag_g'],
        mag_r=tab['mag_r'],
        ra0=ra0,
        dec0=dec0,
        pmra=tab['pmra'],
        pmdec=tab['pmdec'],
        plx=tab['parallax'],
        eplx=tab['parallax_error'],
        epmra=tab['pmra_error'],
        epmdec=tab['pmdec_error'],
        ruwe=tab['ruwe'],
        SelParam=SelParam,
        SelParamSV=SelParamSV,
    )
    red_set = red_sel(**red_kw)

    if nsv != 0:
        sv_blue_set = blue_sel(mag_g=tab['mag_g'],
                               mag_r=tab['mag_r'],
                               ruwe=tab['ruwe'],
                               sv=True,
                               SelParam=SelParam,
                               SelParamSV=SelParamSV)
        sv_red_set = red_sel(**red_kw, sv=True)
        sv_blue_set = get_rand_subset(sv_blue_set, nsv // 2,
                                      rstate=rstate) & (~blue_set)
        sv_red_set = get_rand_subset(sv_red_set, nsv // 2,
                                     rstate=rstate) & (~red_set)

    # the subset of stars that fit our red selection
    red_priority = np.ones(
        len(tab), dtype=float) * priority_giants  # will only be used for red

    nblue = blue_set.sum()
    nred = red_set.sum()

    #  I assume no overlap between the red/blue
    assert ((blue_set & red_set).sum() == 0)

    if strategy == 'redfrac':
        ngoalred = int(SelParam['GIANTS']['goal_red_frac'] * nfibers)
        ngoalblue = nfibers - ngoalred
        if nblue < ngoalblue and nred < ngoalred:
            n_select_blue = nblue
            n_select_red = nred
        elif nblue < ngoalblue and nred >= ngoalred:
            # no enough blue stars for a target fraction
            n_select_blue = nblue
            n_select_red = ngoalred
        elif nblue >= ngoalblue and nred < ngoalred:
            # no enough red stars for a target fraction
            n_select_blue = nfibers - nred
            n_select_red = nred
        elif nblue >= ngoalblue and nred >= ngoalred:
            # enough blue and red targets so we just stick to targets
            n_select_blue = ngoalblue
            n_select_red = ngoalred
        else:
            raise Exception('cant happen')

    blue_set = get_rand_subset(blue_set, n_select_blue, rstate=rstate)

    if not select_in_bins:
        red_set_under = get_rand_subset(red_set, n_select_red, rstate=rstate)
    else:
        red_set_under = select_bin(red_set,
                                   tab['mag_r'],
                                   n_select_red,
                                   maglim=[
                                       SelParam['GIANTS']['min_r_mag'],
                                       SelParam['GIANTS']['max_r_mag']
                                   ],
                                   nbins=SelParam['GIANTS']['red_nbins'],
                                   rstate=rstate)
    del red_set
    # to avoid misuse

    filler_set = filler_sel(mag_g=tab['mag_g'],
                            mag_r=tab['mag_r'],
                            SelParam=SelParam)
    filler_set = filler_set & (~blue_set) & (~red_set_under)

    dn = nfibers - red_set_under.sum() - blue_set.sum()
    if dn >= 0:
        filler_set[rstate.permutation(np.nonzero(filler_set)[0])[dn:]] = False
    else:
        filler_set[:] = False

    retsel = (red_set_under | blue_set | filler_set)
    assert (retsel.sum() <= nfibers)
    ret = dict(blue_set=blue_set,
               red_set=red_set_under,
               filler_set=filler_set,
               red_priority=red_priority)
    if nsv != 0:
        ret['extra_sv_blue_set'] = sv_blue_set
        ret['extra_sv_red_set'] = sv_red_set
    return ret


def column_process(T):
    T.rename_column('objid', 'PS1_ID')
    T.rename_column('source_id', 'SOURCE_ID')
    T['SOURCE_ID'][T['SOURCE_ID'] == -1] = 0  # this is because
    T['PS1_ID'][T['PS1_ID'] == -1] = 0  # of different convention
    T['GAIA_REV_ID'] = 3
    # TODO (should be stored inside input files)


def procone(infits,
            target_density_per_fov=600,
            foot=None,
            mask=None,
            target_extra_sv=0,
            rstate=None,
            extinction_correct=True,
            SelParam=None,
            SelParamSV=None):
    """
    Process one input fits file, select stars with the goal target density
    Return mask
    """
    hdr = fitsio.read_header(infits)
    # pyfits.getheader(infits)
    hpx = hdr['HPX']
    nside = hdr['NSIDE']
    if foot is not None and not foot.overlaps(nside, hpx):
        return

    tab = atpy.Table(pyfits.getdata(infits))
    tab = tab[good_object(tab)]

    if 'qso_flag' in tab.columns:
        tab = tab[~tab['qso_flag']]
    if mask is not None:
        mask_inside = mask.inside(tab['ra'], tab['dec'])
        tab = tab[~mask_inside]
    column_process(tab)
    ra0, dec0 = healpy.pix2ang(nside, hpx, lonlat=True, nest=True)
    tab1 = atpy.Table(tab, copy=True)

    ebv = SFDQuery.query_equ(tab['ra'], tab['dec'])
    # columns used for selection
    ext_g = 3.172 * ebv
    ext_r = 2.271 * ebv
    # coefficients from schalfly 2011
    if not extinction_correct:
        # useful for mocks
        ext_g = 0
        ext_r = 0
    tab1['mag_g'] = tab['gmeanpsfmag'] - ext_g
    tab1['mag_r'] = tab['rmeanpsfmag'] - ext_r

    # we also extinction correct

    pixel_area = (360**2 / np.pi / 12. / nside**2)

    # numbers per healpuxels
    nfibers_goal = int(pixel_area / weave_fov * target_density_per_fov *
                       SelParam['over_subscription'])
    nfibers_poe_goal = math.ceil(pixel_area * SelParam['POE']['maxdens'])
    nfibers_rpm_goal = math.ceil(pixel_area * SelParam['RPM']['maxdens'])
    D = main_sel(tab1,
                 nfibers=nfibers_goal,
                 nsv=target_extra_sv,
                 rstate=rstate,
                 SelParam=SelParam,
                 SelParamSV=SelParamSV)

    poe_sel = gaiapoe_sel(ebv=ebv,
                          phot_g_mean_mag=tab1['phot_g_mean_mag'],
                          phot_bp_mean_mag=tab1['phot_bp_mean_mag'],
                          phot_rp_mean_mag=tab1['phot_rp_mean_mag'],
                          pmra=tab1['pmra'],
                          pmdec=tab1['pmdec'],
                          parallax=tab1['parallax'],
                          parallax_error=tab1['parallax_error'],
                          ruwe=tab1['ruwe'],
                          Ntarget=nfibers_poe_goal,
                          rstate=rstate,
                          SelParam=SelParam,
                          SelParamSV=SelParamSV)
    rpm_sel = gaiarpm_sel(ebv=ebv,
                          phot_g_mean_mag=tab1['phot_g_mean_mag'],
                          phot_bp_mean_mag=tab1['phot_bp_mean_mag'],
                          phot_rp_mean_mag=tab1['phot_rp_mean_mag'],
                          pmra=tab1['pmra'],
                          pmdec=tab1['pmdec'],
                          pmra_error=tab1['pmra_error'],
                          pmdec_error=tab1['pmdec_error'],
                          pmra_pmdec_corr=tab1['pmra_pmdec_corr'],
                          parallax=tab1['parallax'],
                          parallax_error=tab1['parallax_error'],
                          ruwe=tab1['ruwe'],
                          Ntarget=nfibers_rpm_goal,
                          rstate=rstate,
                          SelParam=SelParam,
                          SelParamSV=SelParamSV)
    blue, red, filler, red_priority = D['blue_set'], D['red_set'], D[
        'filler_set'], D['red_priority']

    tab_red = tab[red]
    tab_red['PRIORITY'] = red_priority[red]
    ret = dict(MSTO=tab[blue],
               GIANTS=tab_red,
               FILLER=tab[filler],
               GAIAPOE=tab[poe_sel],
               GAIARPM=tab[rpm_sel])
    if target_extra_sv != 0:
        ret['EXTRASVMSTO'] = tab[D['extra_sv_blue_set']]
        ret['EXTRASVGIANTS'] = tab[D['extra_sv_red_set']]
    return ret


def procall(inprefix,
            oprefix='./',
            version='333333',
            target_density_per_fov=630,
            target_extra_sv=0,
            seed=14343434,
            footprint_file=None,
            extinction_correct=True,
            selection_config_file=None):
    rstate = np.random.default_rng(seed)

    SelParam, SelParamSV = getParams(selection_config_file)
    print('Creating initial target selection files')
    foot = weave_galr.utils.footprint.read_fields(footprint_file)
    # sort the list so we get determinstic behaviour
    fits = sorted(glob.glob(inprefix + '/*/*/*fits'))

    tabs = {
        'MSTO': [],
        'GIANTS': [],
        'FILLER': [],
        'GAIAPOE': [],
        'GAIARPM': []
    }
    if target_extra_sv != 0:
        tabs['EXTRASVMSTO'] = []
        tabs['EXTRASVGIANTS'] = []
    for f in fits:
        ret = procone(f,
                      target_density_per_fov,
                      foot=foot,
                      mask=SkyMask,
                      target_extra_sv=target_extra_sv,
                      rstate=rstate,
                      extinction_correct=extinction_correct,
                      SelParam=SelParam,
                      SelParamSV=SelParamSV)
        if ret is None:
            continue
        for k in ret.keys():
            tabs[k].append(ret[k])

    for k in list(tabs.keys()):
        tabs[k] = atpy.vstack(tabs[k])

    for curt in tabs.values():
        curt.rename_column('ra', 'RA')
        curt.rename_column('dec', 'DEC')

    if not os.path.exists(oprefix):
        os.makedirs(oprefix)
    for name, curt in tabs.items():
        header = pyfits.Header()
        header['CATALOG'] = name
        header['VERSION'] = version
        fill_selconfig_header(header, SelParam, SelParamSV)
        hdu_list = pyfits.HDUList(
            [pyfits.PrimaryHDU(header=header),
             pyfits.BinTableHDU(curt)])
        hdu_list.writeto(f'{oprefix}/{name}_{version}.fits',
                         overwrite=True,
                         checksum=True)
