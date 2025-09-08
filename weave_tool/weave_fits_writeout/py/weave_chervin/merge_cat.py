import re
import hashlib
import datetime
import astropy.io.fits as pyfits
import astropy.table as atpy
import numpy as np
import sqlutilpy as sqlutil
import dustmaps.sfd
import sqlalchemy
import pandas
import healpy
import weave_chervin.utils as util
import weave_chervin.crossmatcher as crossmatcher
import weave_chervin._version as VE
import weave_chervin.utils.footprint as UF

wsdb = util.get_wsdb_host()


def get_gaia_info(GAIA_REV):
    """return dictionary with info about the gaia revision X
    """
    D = {}
    if GAIA_REV == 3:
        D['GAIA_TABLE'] = 'gaia_dr3.gaia_source'
        D['GAIA_EPOCH'] = 2016
    else:
        raise Exception("shouldn't happen")
    return D


def make_sky(foot, density_per_deg, GAIA_REV=None, rstate=None):
    """
    Produce the catalog of skys
    """

    skynside = 64
    # this is the hpx nside within which I produce skys

    assert (skynside < foot.nside)
    # I need to verify that the sky hpx is larger than the footprint pixels

    GAIA_TABLE, GAIA_EPOCH = [
        get_gaia_info(GAIA_REV)[_] for _ in ['GAIA_TABLE', 'GAIA_EPOCH']
    ]

    foot_hpx = np.unique(foot.hpx // (foot.nside // skynside)**2)
    # unique hpxels to consider in the skynside resolution
    area_pix = 360.**2 / np.pi / 12. / skynside**2
    nsky = int(area_pix * density_per_deg)
    # target number of sky fibers
    maxnside = 1 << 29
    # generate hpx grid points at this resolution
    rat = (maxnside // skynside)**2
    ras = []
    decs = []
    hpxs = []
    aperture = 3  # arcseconds to check for presence of objects
    mult = 3
    # how much to oversample take into account that some locations will
    # be rejected
    print('sky gaia')
    for curhp in foot_hpx:
        ids = rstate.integers(curhp * rat, curhp * rat + rat, size=mult * nsky)
        curras, curdecs = healpy.pix2ang(maxnside, ids, nest=True, lonlat=True)
        ras.append(curras)
        decs.append(curdecs)
        hpxs.append(np.zeros(len(curras), dtype=int) + curhp)
    ras, decs, hpxs = [np.concatenate(_) for _ in [ras, decs, hpxs]]

    xras, = crossmatcher.doit(GAIA_TABLE,
                              ras,
                              decs,
                              'tt.ra',
                              rad=aperture,
                              host=wsdb,
                              db='wsdb')
    object_present = np.isfinite(xras)
    ras, decs, hpxs = (ras[~object_present], decs[~object_present],
                       hpxs[~object_present])
    print('sky ps1')
    xras, = crossmatcher.doit('panstarrs_dr1.stackobjectthin',
                              ras,
                              decs,
                              'tt.ra',
                              rad=aperture,
                              host=wsdb,
                              db='wsdb')
    object_present = np.isfinite(xras)
    ras, decs, hpxs = (ras[~object_present], decs[~object_present],
                       hpxs[~object_present])
    ras_all = []
    decs_all = []
    for curhp in foot_hpx:
        i1, i2 = (np.searchsorted(hpxs, curhp, 'left'),
                  np.searchsorted(hpxs, curhp, 'right'))
        if i1 == i2:
            print('WARNING no sky for hpx')
            continue
        if (i2 - i1) <= nsky:
            print('warning not enough sky', i2 - i1, nsky)
        else:
            i2 = i1 + nsky + 1
        ras_all.append(ras[slice(i1, i2)])
        decs_all.append(decs[slice(i1, i2)])
    ras, decs = [np.concatenate(_) for _ in [ras_all, decs_all]]
    xind = foot.contains_pos(ras, decs)
    ras, decs = [_[xind] for _ in [ras, decs]]
    # TODO check the ids of skys
    tab = atpy.Table({
        'GAIA_RA': ras,
        'GAIA_DEC': decs,
        'GAIA_EPOCH': np.zeros(len(ras)) + GAIA_EPOCH,
        'GAIA_DR': np.zeros(len(ras), dtype=int) + GAIA_REV,
        'SOURCE_ID': np.zeros(len(ras), dtype=int),
        'PS1_ID': np.zeros(len(ras), dtype=int),
        'GAIA_PMRA': np.zeros(len(ras)) + 0,
        'GAIA_PMDEC': np.zeros(len(ras)) + 0,
        'GAIA_PMRA_ERR': np.zeros(len(ras)) + 0,
        'GAIA_PMDEC_ERR': np.zeros(len(ras)) + 0,
        'GAIA_PARAL': np.zeros(len(ras)) + 0,
        'GAIA_PARAL_ERR': np.zeros(len(ras)) + 0,
        'OPTCAT_MAG_G': np.zeros(len(ras)) + np.nan,
        'OPTCAT_MAG_R': np.zeros(len(ras)) + np.nan,
        'OPTCAT_MAG_I': np.zeros(len(ras)) + np.nan,
        'MAG_G': np.zeros(len(ras)) + np.nan,
        'MAG_R': np.zeros(len(ras)) + np.nan,
        'MAG_I': np.zeros(len(ras)) + np.nan,
        'MAG_G_ERR': np.zeros(len(ras)) + np.nan,
        'MAG_R_ERR': np.zeros(len(ras)) + np.nan,
        'MAG_I_ERR': np.zeros(len(ras)) + np.nan,
    })
    return tab


def packer(arrs, bits):
    """ pack arrays into bitmask
    i.e
    packer((arr1,arr2,arr3),(3,4,5))
    where arr1 will take the least significant 3 bit,
    arr2 middle 4
    arr3 most significant
    """
    cumb = 0
    ret = np.zeros(len(arrs[0]), dtype=np.int64)
    for cura, curb in zip(arrs, bits):
        assert (np.all(cura < (1 << curb)))
        ret[:] = ret + (cura << cumb)
        cumb += curb
    return ret


def make_target_id(npt, quarter=1, shared=True, pointed=False):
    """
    Make target id based on sequential number of the object,
    survey and what kind of quarter we are

    I contruct targetid like this
    24/3/6 bits
    """
    if shared and not pointed:
        prog = 0
    elif not shared and pointed:
        prog = 1
    else:
        raise Exception('cant be pointed and shared simultaneously')
    ids = np.arange(npt)
    ids = packer((ids, prog, quarter), (24, 3, 6))
    assert (npt < 2**24)
    return ids


def filler_id_coordinates_all(tabs, objtypes, GAIA_REV=None):
    sizes = [len(tabs[_]) for _ in objtypes]
    cumsizes = np.concatenate([[0], np.cumsum(sizes)])
    combined_tab = atpy.vstack([tabs[_] for _ in objtypes],
                               metadata_conflicts='silent')
    # I'm doing here things by survey to speed thing up
    print('Filling ids ps1')
    new_ps1_id = filler_id_coordinates(combined_tab, 'ps1', GAIA_REV=GAIA_REV)
    print('Filling ids gaia')
    new_gaia_id = filler_id_coordinates(combined_tab,
                                        'gaia',
                                        GAIA_REV=GAIA_REV)
    for i, curt in enumerate(objtypes):
        tabs[curt]['PS1_ID'][:] = new_ps1_id[cumsizes[i]:cumsizes[i + 1]]
        tabs[curt]['SOURCE_ID'][:] = new_gaia_id[cumsizes[i]:cumsizes[i + 1]]
        assert (tabs[curt]['PS1_ID'] >= 0).all()
        assert (tabs[curt]['SOURCE_ID'] >= 0).all()


def filler_id_coordinates(tab, survey, GAIA_REV=None):
    """ This checks the IDS in the input catalog
    and fills the PS1 or GAIA ID if needed
    It updates the catalog in place
    """
    ps1_id = tab['PS1_ID']
    gaia_id = tab['SOURCE_ID']
    ra = tab['RA']
    dec = tab['DEC']
    rad = 1

    gaia_ps_source = (gaia_id > 0) & (ps1_id > 0)
    ps_only_source = (ps1_id > 0)
    gaia_only_source = (gaia_id > 0)

    GAIA_TABLE, GAIA_EPOCH = [
        get_gaia_info(GAIA_REV)[_] for _ in ['GAIA_TABLE', 'GAIA_EPOCH']
    ]

    # CHECK that everything has either GAIA or PS1 ID
    assert (((~ps_only_source) & (~gaia_only_source) &
             (gaia_ps_source)).sum() == 0)
    if survey == 'ps1':
        ups1_id = np.unique(ps1_id[ps1_id > 0])
        xps1_id, = sqlutil.local_join(
            '''
            select p.objid  from panstarrs_dr1.meanobject as p,
            mytab as m where p.objid=m.objid;
                        ''',
            'mytab',
            (ups1_id, ),
            ('objid', ),
            host=wsdb,
            db='wsdb',
        )
        yps1_id, = sqlutil.local_join(
            '''
            select p.objid  from panstarrs_dr1.objectthin as p,  mytab as m
            where p.objid=m.objid;
                        ''',
            'mytab',
            (ups1_id, ),
            ('objid', ),
            host=wsdb,
            db='wsdb',
        )

        # here we check that the ps1 id actually exist
        assert np.all(np.in1d(ups1_id, xps1_id) | np.in1d(ups1_id, yps1_id))

        # Now we try to check if the targets have missing PS1 or Gaia info
        xps1_id, xdist = crossmatcher.doit(
            'panstarrs_dr1.objectthin',
            ra,
            dec,
            'objid, 3600*q3c_dist(tt.ramean,tt.decmean,m.ra,m.dec)',
            host=wsdb,
            db='wsdb',
            radeccols=('ramean', 'decmean'),
            rad=rad)
        assert (np.all(xdist[np.isfinite(xdist)] < 1))

        ps1_id[(ps1_id == 0) & (xps1_id > 0)] = xps1_id[(ps1_id == 0)
                                                        & (xps1_id > 0)]
        return ps1_id
    elif survey == 'gaia':
        xgaia_id, = sqlutil.local_join(f'''
    select s.source_id from mytab m left join
                  {GAIA_TABLE} as s
    on (s.source_id=m.source_id)
                    order by xid;''',
                                       'mytab',
                                       (np.arange(len(gaia_id)), gaia_id),
                                       ('xid', 'source_id'),
                                       host=wsdb,
                                       db='wsdb')
        assert (np.all((xgaia_id == gaia_id)[gaia_id > 0]))

        xgaia_id, xdist = crossmatcher.doit(
            GAIA_TABLE,
            ra,
            dec,
            'source_id, 3600*q3c_dist(tt.ra,tt.dec,m.ra,m.dec)',
            host=wsdb,
            db='wsdb',
            rad=rad)
        assert (np.all(xdist[np.isfinite(xdist)] < 1))
        gaia_id[gaia_id > 0] = xgaia_id[gaia_id > 0]
        return gaia_id
    else:
        raise Exception('oops')


def fetch_gaia_info(gaia_id, GAIA_REV=None):
    """
    Get various columns from gaia
    """
    ret = {}
    GAIA_TABLE, GAIA_EPOCH = [
        get_gaia_info(GAIA_REV)[_] for _ in ['GAIA_TABLE', 'GAIA_EPOCH']
    ]

    D = sqlutil.local_join(f'''
select pmra,pmdec,pmra_error,pmdec_error,
parallax, parallax_error,
ra,dec,ref_epoch,
phot_bp_mean_mag,
phot_rp_mean_mag,
phot_g_mean_mag,
phot_bp_mean_flux_over_error,
phot_rp_mean_flux_over_error,
phot_g_mean_flux_over_error,
l, b
  from mytab m
left join
              {GAIA_TABLE} as s
on (s.source_id=m.source_id)
                order by xid;''',
                           'mytab', (np.arange(len(gaia_id)), gaia_id),
                           ('xid', 'source_id'),
                           host=wsdb,
                           db='wsdb',
                           asDict=True)
    ret['GAIA_RA'] = D['ra']
    ret['GAIA_DEC'] = D['dec']
    ret['GAIA_EPOCH'] = D['ref_epoch']
    ret['GAIA_PMRA'] = D['pmra']
    ret['GAIA_PMDEC'] = D['pmdec']
    ret['GAIA_PMRA_ERR'] = D['pmra_error']
    ret['GAIA_PMDEC_ERR'] = D['pmdec_error']
    ret['GAIA_PARAL'] = D['parallax']
    ret['GAIA_PARAL_ERR'] = D['parallax_error']
    ret['GAIA_MAG_G'] = D['phot_g_mean_mag']
    ret['GAIA_MAG_BP'] = D['phot_bp_mean_mag']
    ret['GAIA_MAG_RP'] = D['phot_rp_mean_mag']

    for filt in ['g', 'bp', 'rp']:
        curerr = 2.5 / np.log(10) / D['phot_' + filt + '_mean_flux_over_error']
        curerr[curerr > 1] = np.nan
        ret['GAIA_MAG_' + filt.upper() + '_ERR'] = curerr
    return ret


def fetch_ps1_phot(objid):
    """
    fetch ps1 photometry
    """
    ret = {}
    D = sqlutil.local_join('''
SELECT x.* from mytab AS m LEFT JOIN LATERAL (
    SELECT ramean,decmean,
        gmeanpsfmag,
        gmeanpsfmagerr,
        rmeanpsfmag,
        rmeanpsfmagerr,
        imeanpsfmag,
        imeanpsfmagerr
        FROM
            panstarrs_dr1.objectthin AS so
        LEFT JOIN
            panstarrs_dr1.meanobject AS s ON ( so.objid=s.objid)
        WHERE so.objid=m.objid
        LIMIT 1
   ) AS x ON true
        ORDER BY xid;''',
                           'mytab', (np.arange(len(objid)), objid),
                           ('xid', 'objid'),
                           host=wsdb,
                           db='wsdb',
                           asDict=True)
    assert (len(D['gmeanpsfmag']) == len(objid))
    for col in 'gri':
        colU = col.upper()
        colval = D[f'{col}meanpsfmag']
        ecolval = D[f'{col}meanpsfmagerr']
        xind = (ecolval < 0) | (ecolval > 1)
        colval[xind] = np.nan
        ecolval[xind] = np.nan
        ret[f'PS_MAG_{colU}'] = colval
        ret[f'PS_MAG_{colU}_ERR'] = ecolval
        ret[f'MAG_{colU}'] = colval
        ret[f'MAGERR_{colU}'] = ecolval
        ret[f'OPTCAT_MAG_{colU}'] = colval
        ret[f'OPTCAT_MAGERR_{colU}'] = ecolval
    ret['PS_DR'] = 1
    ret['OPTCAT'] = 'PanStarrs'
    ret['OPTCAT_DR'] = 1
    ret['PS_ID'] = objid
    return ret, D['ramean'], D['decmean']


def fetch_ebv(ra, dec):
    """
    fetch sfd ebv
    """
    ret = {}
    Q = dustmaps.sfd.SFDQuery()
    ebv = Q.query_equ(ra, dec)
    ret['EBVCAT'] = 'Schlegel'
    ret['EBVCAT_DR'] = '1998'
    ret['EBVCAT_RA'] = ra
    ret['EBVCAT_DEC'] = dec
    ret['EBVCAT_EBV'] = ebv
    ret['EBVCAT_EBV_ERR'] = np.nan
    return ret


def fetch_sdss_phot(ra, dec):
    """
    fetch sdss magnitudes
    """
    ret = {}
    D = sqlutil.local_join('''
select x.*
from mytab m
left join lateral( select s.ra,s.dec, objid,
psfmag_u, psfmagerr_u,
psfmag_g, psfmagerr_g,
psfmag_r, psfmagerr_r,
psfmag_i, psfmagerr_i,
psfmag_z, psfmagerr_z from
              sdssdr14.photoobjall as s
where q3c_join(m.ra,m.dec,s.ra,s.dec,1./3600)
order by (mode=1, type=6) desc limit 1)
as x on true
                order by xid;''',
                           'mytab', (np.arange(len(ra)), ra, dec),
                           ('xid', 'ra', 'dec'),
                           host=wsdb,
                           db='wsdb',
                           asDict=True)
    assert (len(ra) == len(D['objid']))
    for col in 'ugriz':
        colval = D[f'psfmag_{col}']
        ecolval = D[f'psfmagerr_{col}']
        colU = col.upper()
        xind = (ecolval < 0) | (ecolval > 1) | (colval < 0)
        colval[xind] = np.nan
        ecolval[xind] = np.nan
        ret[f'SDSS_MAG_{colU}'] = colval
        ret[f'SDSS_MAG_{colU}_ERR'] = ecolval
    ret['SDSS_DR'] = 14
    ret['SDSS_ID'] = D['objid']
    ret['SDSS_RA'] = D['ra']
    ret['SDSS_DEC'] = D['dec']
    return ret


def convert_priority(x, mi, ma):
    """
    rescale priorities from 1-10 to mi,ma range
    """
    return mi + ((x - 1) / 9) * (ma - mi)


def convert_tab(tab, objtype, GAIA_REV=None):
    """
    This takes the input catalog, renames the relevant columns
if needed, scales the priorities and drops the not used columns
    """
    mapperD = {
        'PS1_ID': 'PS_ID',
        'SOURCE_ID': 'GAIA_ID',
        'GAIA_REV_ID': 'GAIA_DR',
        'RA': 'RA_0',
        'DEC': 'DEC_0'
    }
    copy_cols = [
        'FEHCAT', 'FEHCAT_DEC', 'FEHCAT_DR', 'FEHCAT_FEH', 'FEHCAT_FEH_ERR',
        'FEHCAT_ID', 'FEHCAT_RA'
    ]
    ret = {}
    for k in mapperD.keys():
        curcol = tab[k]

        if k == 'GAIA_REV_ID':
            curcol = curcol.astype(np.int64)
            assert (np.all((curcol == GAIA_REV) | (curcol == 0)))
        if k == 'PS1_ID' or k == 'SOURCE_ID':
            # because the validator tests it
            curcol = curcol.astype(np.int64)
        if k == 'RA' or k == 'DEC':
            # because the validator tests it
            curcol = curcol.astype(np.float64)
        ret[mapperD[k]] = curcol
    mi, ma = 5, 5
    for k in copy_cols:
        if k in tab.columns:
            ret[k] = tab[k]
    if 'PRIORITY' in tab.columns:
        ret['TARGPRIO'] = tab['PRIORITY']#convert_priority(tab['PRIORITY'].astype(np.float64), mi, ma)
    else:
        assert (mi == ma)
        ret['TARGPRIO'] = np.zeros(len(tab), dtype=np.float64) + mi
    return atpy.Table(ret)


def get_obstemp(seeing=1,
                transparency=.8,
                airmass=1.5,
                moon=0,
                brightness=21.5):
    """
    convert observational conditions into obstemp string
    """
    x1 = {7: 'A', 8: 'B', 9: 'C', 10: 'D'}[int(seeing * 10)]
    x2 = {
        9: 'A',
        8: 'A',
        7: 'B',
        6: 'C',
        5: 'D',
        4: 'E'
    }[int(transparency * 10)]
    x3 = {13: 'A', 14: 'B', 15: 'C', 16: 'D'}[int(airmass * 10)]
    x4 = {90: 'A', 70: 'B', 50: 'C', 30: 'D', 0: 'E'}[moon]
    x5 = {217: 'A', 215: 'B', 210: 'C'}[int(brightness * 10)]
    return (x1 + x2 + x3 + x4 + x5)


def bitmask2targprog(bitmask, bitmD, nameD):
    """
    Convert the bitmask integer to progtemp like string
    """
    res = []
    delim = '|'
    for k, v in bitmD.items():
        if (bitmask & v) > 0:
            res.append(nameD[k])
    return delim + delim.join(res) + delim


def bitmask2aps(bitmask, bitmD, apsD):
    """
    Convert the bitmask integer to aps like string
    """
    curl = [apsD['__STAR__']]  # always use APS STAR flags
    for k, v in bitmD.items():
        if (bitmask & v) > 0:
            curl.append(apsD[k])
    arr = np.array(curl).sum(axis=0) > 0
    return ''.join((arr > 0).astype(int).astype(str))


def get_targprog(tab, shared=True, pointed=False):
    """
    Compute TARGPROG strings
    """
    print('getting TARGPROG')
    maxLen = 40
    res = []
    bitmD = util.get_object_bitmask_dict()
    nameD = util.get_object_targprog_dict()
    if shared:
        survey = 'SHA|'
    elif pointed:
        survey = 'POI|'
    else:
        raise Exception("shouldn't happen")
    for i in range(len(targbits)):
        res.append(bitmask2targprog(targbits[i], bitmD, nameD) + survey)
        assert (len(res[-1]) < maxLen)
    return np.array(res)


def get_aps_flag(tab):
    """
    Compute the APS flags
    """
    print('getting APS')
    targbits = tab['GA_TARGBITS']
    bitmD = util.get_object_bitmask_dict()
    apsD = util.get_aps_flag_dict()
    res = []
    for i in range(len(targbits)):
        res.append(bitmask2aps(targbits[i], bitmD, apsD))
    return np.array(res)


def md5_sum(fname):
    m = hashlib.md5()
    with open(fname, 'rb') as fp:
        m.update(fp.read())
    return m.hexdigest()


def validate_trimester(trimester):
    """ convert trimester 2021A1 to integer number """
    M = re.match('(20[0-9]{2})([AB])([12])', trimester)
    if M is None:
        raise Exception('incorrect trimester format')
    year = int(M.group(1))
    assert (year >= 2021)
    trim_letter = M.group(2)
    trim_number = int(M.group(3))
    quarter = (year - 2022) * 4 + {
        'A': 0,
        'B': 1
    }[trim_letter] * 2 + (trim_number - 1)
    assert (quarter >= 0)
    return quarter


def doit(cats,
         progtemps=None,
         program=None,
         footprint_file=None,
         targcat=None,
         targsurvey=None,
         GAIA_REV=None,
         trimester=None):
    """
    Main program that creates the combined file
    """
    if program == 'shared':
        shared = True
    elif program == 'pointed':
        shared = False
    else:
        raise Exception("can't happen")
    pointed = not shared

    foot = None
    if footprint_file is not None:
        foot = UF.read_fields(footprint_file)

    GAIA_TABLE, GAIA_EPOCH = [
        get_gaia_info(GAIA_REV)[_] for _ in ['GAIA_TABLE', 'GAIA_EPOCH']
    ]

    if targsurvey is None:
        targsurvey = 'GA-LRHIGHLAT'

    # https://docs.google.com/document/d/1NDFni1iVHtr4GOA5Uvc-RtdvzusNbW_bpw4xcEufegE/edit
    obstemp = 'DACEB'
    quarter = validate_trimester(trimester)

    sky_density = util.get_sky_density()  # per degree
    if targcat is None:
        targcat = f'GA-LRHIGHLAT_{trimester}.fits'

    meta = {}
    tabs = {}
    sizes = {}
    coltypes = {}
    counter = {}
    cata_meta = []
    for ii, curcat in enumerate(cats):
        tab = atpy.Table().read(curcat)
        hdr = pyfits.getheader(curcat, 1)
        hdr0 = pyfits.getheader(curcat)
        objtype = hdr.get('CATALOG') or hdr0.get('CATALOG')
        print('Processing', objtype, curcat)
        if len(tab) == 0:
            print("WARNING %s is empty" % (objtype))
        assert (tab['SOURCE_ID'].min() >= 0)
        cata_meta.append({
            'file': curcat,
            'objtype': objtype,
            'md5': md5_sum(curcat),
        })
        for k in hdr0.keys():
            if re.match('.*Selectio.*', k) is not None:
                if k in meta:
                    assert hdr0[k] == meta[k]
                else:
                    meta[k] = hdr0[k]

        if foot is not None:
            xind = foot.contains_pos(tab['RA'], tab['DEC'])
            tab = tab[xind]
        hpx = healpy.ang2pix(1 << 29,
                             tab['RA'],
                             tab['DEC'],
                             nest=True,
                             lonlat=True)
        tab = tab[np.argsort(hpx)]
        objtype = 'DUMMY'
        if objtype in counter:
            raise Exception(
                " The file list has more than one catalog with object type " +
                objtype)
            # ensure that we don't have more than one catalog of
            # a given type
        counter[objtype] = 1
        tabs[objtype] = tab
        sizes[objtype] = len(tab)
        assert np.issubdtype(tab['SOURCE_ID'].dtype, np.int64)
        assert np.issubdtype(tab['PS1_ID'].dtype, np.int64)
        del tab
    objtypes = sorted(list(tabs.keys()), key=lambda x: -sizes[x])
    filler_id_coordinates_all(tabs, objtypes, GAIA_REV=GAIA_REV)

    for objtype in objtypes:
        tab = convert_tab(tabs[objtype], objtype, GAIA_REV=GAIA_REV)
        tabs[objtype] = tab
        for k in tab.columns:
            if k not in coltypes:
                coltypes[k] = tab[k].dtype
            try:
                assert (coltypes[k] == tab[k].dtype)
            except:  # noqa
                print('type mismatch', k, coltypes[k], tab[k].dtype)
                raise
        del tab
    # now we start merging

    print('merging')
    # start from the largest
    curcat = tabs[objtypes[0]].to_pandas()
    curcat['rowid'] = np.arange(len(curcat))
    maxid = len(curcat)
    eng = sqlalchemy.create_engine('sqlite:///')
    curcat.to_sql('curcat', eng)
    with eng.connect() as conn:
        conn.execute(
            sqlalchemy.text(
                'update curcat set "PS_ID" = null where "PS_ID"=0'))
        conn.execute(
            sqlalchemy.text(
                'update curcat set GAIA_ID = null where GAIA_ID=0'))
        conn.execute(sqlalchemy.text('create index ii1 on curcat(GAIA_ID)'))
        conn.execute(sqlalchemy.text('create index ii2 on curcat("PS_ID")'))
        conn.execute(sqlalchemy.text('analyze'))
    # these are the current columns we have in curcat
    curcolumns = [_ for _ in curcat.columns]

    # here we iterate over all the catalogs
    # in each iteration we update the matching records
    # to coadd targetbits and fill unfilled columns
    # and then we append the new records
    for curtype in objtypes[1:]:
        print('merging', curtype)
        newt = tabs[curtype].to_pandas()
        # add primary key
        newt['rowid'] = maxid + np.arange(len(newt))
        maxid = newt['rowid'].max()
        # upload the data
        newt.to_sql('newcat', eng)
        # replace zeros by nulls
        eng.execute('update newcat set "PS_ID" = null where "PS_ID"=0')
        eng.execute('update newcat set GAIA_ID = null where GAIA_ID=0')
        eng.execute('create index ii5 on newcat(GAIA_ID)')
        eng.execute('create index ii6 on newcat("PS_ID")')
        eng.execute('create index ii7 on newcat(rowid)')
        eng.execute('analyze')
        # here we create the match table to identify  the targets
        # that already exist

        # stars with matching gaia id and no ps id
        eng.execute('''create table matches as select c.rowid as rowid1 ,
n.rowid as rowid2 from curcat as c, newcat as n
where c."PS_ID" is null and
 n.PS_ID is null
and c.GAIA_ID = n.GAIA_ID
''')
        # stars with matching ps id and no gaia id
        eng.execute('''insert into matches select c.rowid as rowid1 ,
n.rowid as rowid2 from curcat as c, newcat as n
where c."GAIA_ID" is null and
 n.GAIA_ID is null
and c.PS_ID = n.PS_ID
''')
        # stars with matching ps id and gaia id
        eng.execute('''insert into matches select c.rowid as rowid1 ,
n.rowid as rowid2 from curcat as c, newcat as n
where c."GAIA_ID"=n.GAIA_ID
and c.PS_ID = n.PS_ID
''')
        eng.execute('create index ii3 on matches(rowid1)')
        eng.execute('create index ii4 on matches(rowid2)')
        eng.execute('analyze')

        # these are the columns that we have to add
        newcolumns = [_ for _ in newt.columns if _ not in curcolumns]

        # the net table with have old columns + new columns
        collist = ','.join(['c.' + _ for _ in curcolumns] +
                           ['n.' + _ for _ in newcolumns])

        eng.execute(
            f'''create table combined as select {collist} from newcat as n , curcat as c
 limit 0
''')

        # copy existing data
        collist = ','.join([_ for _ in curcolumns])
        eng.execute(
            f'''insert into combined ({collist}) select {collist} from curcat as c
''')

        # update the columns from the new table for existing targets
        # do not *YET* update targeting bits and priorities

        for curcol in curcolumns:
            if curcol in [
                    'PS_ID', 'GAIA_ID', 'rowid', 'GA_TARGBITS', 'TARGPRIO'
            ]:
                continue
            if curcol in newt.columns:
                eng.execute(f'''
UPDATE combined SET {curcol} =
    (
    SELECT coalesce(n.{curcol},combined.{curcol})
        from newcat as n,
             matches as m
        where n.rowid=m.rowid2 and m.rowid1=combined.rowid
    )
    where exists (select 1 from matches as m where m.rowid1=rowid)
''')
        # combine targeting bits
        eng.execute('''update combined set GA_TARGBITS = (
select combined.GA_TARGBITS +  n.GA_TARGBITS from newcat as n, matches as m
where n.rowid = m.rowid2  and m.rowid1=combined.rowid
) where exists ( select 1 from matches as m where m.rowid1 = combined.rowid)
''')
        # combine targeting priorities
        eng.execute('''update combined set TARGPRIO = (
        select max(combined.TARGPRIO, n.TARGPRIO) from newcat as n,
matches as m
where n.rowid = m.rowid2  and m.rowid1=combined.rowid
) where exists ( select 1 from matches as m where m.rowid1 = combined.rowid)
''')
        collist = ','.join([_ for _ in newt.columns])
        eng.execute(
            f'''insert into combined({collist}) select {collist} from newcat
where not exists(select 1 from matches as m where m.rowid2=newcat.rowid)''')
        eng.execute('drop table curcat')
        eng.execute('create table curcat as select * from combined')
        eng.execute('create index ii1 on curcat(GAIA_ID)')
        eng.execute('create index ii2 on curcat("PS_ID")')
        eng.execute('analyze')
        eng.execute('drop table combined')
        eng.execute('drop table matches')
        eng.execute('drop table newcat')
        curcolumns = set(curcolumns).union(set(newcolumns))
    collist = []
    maskv = -9999

    for curc in curcolumns:
        if curc == 'rowid':
            continue
        if np.issubdtype(coltypes[curc], int):
            collist.append(f'coalesce({curc},{maskv}) as {curc}')
        else:
            collist.append(curc)
    collist = ','.join(collist)
    pdTab = pandas.read_sql(f'select {collist} from curcat',
                            eng,
                            coerce_float=False)
    resTab = {}
    for k in pdTab.columns:
        curarr = pdTab[k]
        if np.issubdtype(coltypes[k], int):
            resTab[k] = np.ma.array(curarr, mask=curarr == maskv)
        if np.issubdtype(coltypes[k], np.floating):
            resTab[k] = np.ma.array(
                curarr,
                mask=~np.isfinite(np.array(curarr, dtype=np.float64)),
                dtype=coltypes[k])
        elif np.issubdtype(coltypes[k], np.bytes_):
            curarr = np.array(curarr.str.decode('ascii'), dtype=str)
            resTab[k] = np.ma.array(curarr, mask=curarr == 'None')
        else:
            resTab[k] = curarr

    tab = atpy.Table(resTab)
    # sort by position
    assert (tab['RA_0'].mask.sum() == 0)
    assert (tab['DEC_0'].mask.sum() == 0)
    tab['RA_0'] = (tab['RA_0'] + 360) % 360
    # fix if someone provided broken ras  - naughty-naughty
    hpx = healpy.ang2pix(1 << 29,
                         tab['RA_0'].filled(),
                         tab['DEC_0'].filled(),
                         lonlat=True,
                         nest=True)

    tab = tab[np.argsort(hpx)]
    # I kept those here just for sorting purposes
    # I'd be good to verify that the coordinates match the
    # newly minted GAIA_RA
    tab.remove_columns(['RA_0', 'DEC_0'])

    print('fetching gaia')
    ret = fetch_gaia_info(tab['GAIA_ID'], GAIA_REV=GAIA_REV)
    for k in ret.keys():
        tab[k] = ret[k]

    print('fetching ps1')
    ret, psra, psdec = fetch_ps1_phot(tab['PS_ID'])
    for k in ret.keys():
        tab[k] = ret[k]
    xind = ~np.isfinite(tab['GAIA_RA'])

    # fill coordinates with Gaia info
    tab['GAIA_EPOCH'][xind] = GAIA_EPOCH
    tab['GAIA_DR'][xind] = GAIA_REV
    tab['GAIA_RA'][xind] = psra[xind]
    tab['GAIA_DEC'][xind] = psdec[xind]
    tab['GAIA_PMRA'][xind] = np.nan
    tab['GAIA_PMDEC'][xind] = np.nan
    tab['GAIA_PMRA_ERR'][xind] = np.nan
    tab['GAIA_PMDEC_ERR'][xind] = np.nan
    tab['GAIA_PARAL'][xind] = np.nan
    tab['GAIA_PARAL_ERR'][xind] = np.nan

    ra, dec = tab['GAIA_RA'], tab['GAIA_DEC']
    assert (np.isfinite(np.sum(ra)))
    print('fetching sdss')
    ret = fetch_sdss_phot(ra, dec)
    for k in ret.keys():
        tab[k] = ret[k]
    ret = fetch_ebv(ra, dec)
    for k in ret.keys():
        tab[k] = ret[k]
    nrows = len(tab)
    tab['TARGUSE'] = np.repeat('T', nrows)
    tab['TARGCLASS'] = np.repeat('STAR', nrows)
    # tab['TARGPROG'] = get_targprog(tab, shared=shared, pointed=pointed)
    # tab['APS_FLAG'] = get_aps_flag(tab)

    print('making sky')
    # ensure skies are allocated differently in different quarters
    #rstate = np.random.default_rng(4336360 + quarter)
    # skytab = make_sky(foot, sky_density, GAIA_REV=GAIA_REV, rstate=rstate)
    #sky_prio = util.get_priority('SKY')[0]
    #skytab['TARGPRIO'] = np.zeros(len(skytab), dtype=float) + sky_prio
    #skytab['TARGUSE'] = np.repeat('S', len(skytab))
    # skytab['TARGCLASS'] = np.repeat('SKY', len(skytab))

    #tab = atpy.vstack((tab, skytab))
    #tab['GAIA_ID'] = np.ma.array(tab['GAIA_ID'], mask=tab['GAIA_ID'] <= 0)
    nrows = len(tab)

    tab['TARGSRVY'] = np.repeat(targsurvey, nrows)
    targid = make_target_id(len(tab),
                            quarter=quarter,
                            shared=shared,
                            pointed=pointed)
    tab['TARGID'] = targid
    tab['TARGCAT'] = targcat
    tab['TARGNAME'] = np.array([''] * nrows)
    tab['OBSTEMP'] = np.repeat(obstemp, nrows)
    date = datetime.datetime.now().isoformat()
    ver = VE.version
    # tab["GA_TARGREV"] = ver
    # tab["GA_TARGDATE"] = date
    #if tacalloc is not None:
    #    meta['TACALLOC'] = tacalloc
    bigtab = []
    for curprogtemp in progtemps:
        tab['PROGTEMP'] = np.repeat(curprogtemp, nrows)
        bigtab.append(tab.copy())
    bigtab = atpy.vstack(bigtab)
    return bigtab, meta, cata_meta
