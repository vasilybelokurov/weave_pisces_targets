import os
import astropy.io.fits as pyfits
import numpy as np


def getucds(head, coldefs):
    ret = {}
    for i, c in enumerate(coldefs):
        k = 'TUCD%d' % (i + 1)
        if k in head:
            ret[c.name] = head[k]
    return ret


def getprop(head, coldefs):
    ret = {}
    for i, c in enumerate(coldefs):
        k = 'TPROP%d' % (i + 1)
        if k in head:
            ret[c.name] = head[k]
    return ret


def write_conv_file(tab0,
                    meta,
                    cata_meta,
                    ofile,
                    trimester=None,
                    targsurvey=None):
    dirpref = os.path.dirname(__file__)
    hdu0, hdu1 = pyfits.open(
        f'{dirpref}/data/{targsurvey}_CatalogueTemplate.fits')
    ucds_ref = getucds(hdu1.header, hdu1.columns)
    column_list = hdu1.data.columns

    assert (trimester is not None)
    R1 = []
    nodata_cols = []
    for k in hdu1.columns:
        curdisp = k.disp
        if k.name in tab0.columns:
            if k.format == 'A':
                if k.name == 'PROGTEMP':
                    curform = '7A'
                    curdisp = 'A7'
                else:
                    curform = k.disp
                    curdisp = k.disp
            else:
                curform = k.format
                curdisp = k.disp
            # print(k.name, k.disp, k.format)
            curd = tab0[k.name]
            if isinstance(curd, np.ma.masked_array):
                if k.null is not None:
                    curd = curd.filled(k.null)
                else:
                    if np.issubdtype(curd.dtype, np.floating):
                        curd = curd.filled(np.nan)
                    if k.format[-1] == 'A':
                        curd = curd.astype(str).filled('')
            curk = pyfits.Column(array=curd,
                                 format=curform,
                                 disp=curdisp,
                                 unit=k.unit,
                                 name=k.name,
                                 null=k.null)
        else:
            curk = k
            if k.null is not None:
                nodata_cols.append((k.name, k.null))
            elif k.format[0] in ['E', 'D']:
                nodata_cols.append((k.name, np.nan))
        R1.append(curk)
    # fill columns that I don't provide with null values
    hdu1x = pyfits.BinTableHDU.from_columns(R1)
    for k, v in nodata_cols:
        hdu1x.data[k][:] = v

    for i, c in enumerate(hdu1.columns):
        hdu1x.header['TUCD%d' % (i + 1)] = ucds_ref.get(c.name)
        k = 'TLMIN%d' % (i + 1)
        if k in hdu1.header:
            hdu1x.header[k] = hdu1.header[k]
        k = 'TLMAX%d' % (i + 1)
        if k in hdu1.header:
            hdu1x.header[k] = hdu1.header[k]

        hdu1x.header['TPROP%d' % (i + 1)] = 0
    hdu1x.header['EXTNAME'] = 'WS2025B2-028 CATALOGUE'

    hdu0.header["DATAMVER"] = '8.00'
    hdu0.header['TRIMESTE'] = trimester
    hdu0.header['MAG_G_CM'] = 'OPTCAT_MAG_G'
    hdu0.header['MAG_R_CM'] = 'OPTCAT_MAG_R'
    hdu0.header['MAG_I_CM'] = 'OPTCAT_MAG_I'
    hdu0.header['CAT_MAIL'] = 'alis.j.deason@durham.ac.uk'
    hdu0.header['CAT_CC'] = 'skoposov@ed.ac.uk'
    hdu0.header['CAT_NME1'] = "Alis"
    hdu0.header['CAT_NME2'] = "Deason"
    #TRIMESTE= '2023B2  '           / Observing Trimester                           #TACALLOC= 'WS2023B2I0378'      / Proposal identifier                           TACID   = 'ITP2023-07'         / TAC identifiers                               # ensure the order is the same as in the template
    # Idiotic requirement of WASP if you ask me
    new_header = hdu1.header
    for k, v in hdu1x.header.items():
        new_header[k] = v
    hdu1x.header = new_header
    for k, v in meta.items():
        if len(k) > 8:
            hdu1x.header['HIERARCH ' + k] = v
        else:
            hdu1x.header[k] = v
    hdu1x.header['LONGSTRN'] = 'OGIP 1.0'
    hdul = pyfits.HDUList([hdu0, hdu1x])
    hdul.writeto(ofile, overwrite=True, checksum=True)
