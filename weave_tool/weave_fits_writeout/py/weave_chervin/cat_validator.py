import astropy.io.fits as pyfits
import astropy.table as atpy
import numpy as np
import weave_galr.utils as util

types = util.get_object_types()

allowed_columns = [
    'RA',
    'DEC',
    'SOURCE_ID',
    'PS1_ID',
    'GAIA_REV_ID',
    # metallicies for PRISTINE
    'FEHCAT',
    'FEHCAT_DEC',
    'FEHCAT_DR',
    'FEHCAT_FEH',
    'FEHCAT_FEH_ERR',
    'FEHCAT_ID',
    'FEHCAT_RA',
    'PRIORITY'
]


def validate(fname):
    """ Check the input catalog
    """
    hdr1 = pyfits.getheader(fname, 1)
    hdr = pyfits.getheader(fname)
    # check header
    cat = hdr1.get('CATALOG') or hdr.get('CATALOG')
    version = hdr1.get('VERSION') or hdr.get('VERSION')
    assert (cat is not None)
    assert (version is not None)
    assert (cat in types)
    tab = atpy.Table().read(fname)

    # Check that correct columns exist
    assert ('SOURCE_ID' in tab.columns)
    assert ('RA' in tab.columns)
    assert ('DEC' in tab.columns)
    assert ('GAIA_REV_ID' in tab.columns)
    assert ('PS1_ID' in tab.columns)
    
    # check types
    assert (np.issubdtype(tab['PS1_ID'].dtype, np.int64))
    assert (np.issubdtype(tab['SOURCE_ID'].dtype, np.int64))
    if 'PRIORITY' in tab.columns:
        assert (np.issubdtype(tab['PRIORITY'].dtype, np.floating)
                or np.issubdtype(tab['PRIORITY'].dtype, np.int))

    assert np.all(np.in1d(tab['GAIA_REV_ID'], [0, 3]))
    # gaia dr2 or dr3 or neither, but not mix/match

    # check that source_ids are not too small
    if (tab['SOURCE_ID'] > 0).sum() > 0:
        assert (np.min(tab['SOURCE_ID'][tab["SOURCE_ID"] > 0]) > 4295806720)

    # no negatives
    assert (np.all(tab['PS1_ID'] >= 0))
    assert (np.all(tab['SOURCE_ID'] >= 0))

    # either gaia or ps1 id set
    assert (np.all((tab['SOURCE_ID'] > 0) | (tab['PS1_ID'] > 0)))

    # check ra, dec are goos
    assert (np.all(np.isfinite(tab['RA'])))
    assert (np.all(np.isfinite(tab['DEC'])))
    assert ((tab['RA'] > 0).all())
    assert ((tab['RA'] < 360).all())
    assert ((tab['DEC'] > -90).all())
    assert ((tab['DEC'] < 90).all())

    # check uniqueness
    good_id = tab['SOURCE_ID'][tab['SOURCE_ID'] > 0]
    goodps1_id = tab['PS1_ID'][tab['PS1_ID'] > 0]
    assert (len(good_id) == len(np.unique(good_id)))
    assert (len(goodps1_id) == len(np.unique(goodps1_id)))
    assert (all([_ in allowed_columns for _ in tab.columns]))
