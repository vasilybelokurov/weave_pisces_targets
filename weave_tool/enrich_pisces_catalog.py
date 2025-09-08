#!/usr/bin/env python3
"""
Enrich Pisces target catalog with wsdb data for WEAVE compatibility.

This script takes the basic Pisces catalog (SOURCE_ID, RA, DEC, PRIORITY) 
and enriches it with Gaia EDR3 photometry, astrometry, and cross-match data
required by the WEAVE chervin tool.
"""

import numpy as np
import astropy.io.fits as fits
from astropy.table import Table, join
import sqlutilpy
import os
from tqdm import tqdm

def query_gaia_data(source_ids, batch_size=1000):
    """Query Gaia EDR3 data for given source IDs from wsdb."""
    
    # Split into batches to avoid query size limits
    n_batches = (len(source_ids) + batch_size - 1) // batch_size
    all_results = []
    
    print(f"Querying Gaia data in {n_batches} batches of {batch_size}")
    
    for i in tqdm(range(n_batches), desc="Querying wsdb"):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, len(source_ids))
        batch_ids = source_ids[start_idx:end_idx]
        
        # Convert to comma-separated string for SQL IN clause
        id_list = ','.join(map(str, batch_ids))
        
        query = f"""
        SELECT 
            source_id,
            ra, dec,
            phot_g_mean_mag,
            phot_bp_mean_mag, 
            phot_rp_mean_mag,
            pmra, 
            pmdec, 
            parallax, 
            l as gal_long, 
            b as gal_lat
        FROM gaia_edr3.gaia_source 
        WHERE source_id IN ({id_list})
        """
        
        try:
            batch_result = sqlutilpy.get(query, asDict=True)  # Get as dictionary
            if len(batch_result) > 0:
                all_results.append(batch_result)
        except Exception as e:
            print(f"Error in batch {i}: {e}")
            continue
    
    if all_results:
        # Concatenate all batch results
        import pandas as pd
        df_list = [pd.DataFrame(result) for result in all_results]
        combined_df = pd.concat(df_list, ignore_index=True)
        return combined_df
    else:
        return None

def map_priority_to_target_type(priority):
    """Map Pisces priority to WEAVE target type."""
    mapping = {
        10: 'GIANTS',    # Clean RGB Giants
        9:  'GIANTS',    # Additional RGB Giants  
        8:  'BHB',       # BHB + RR Lyrae (needs sub-classification)
        7:  'GIANTS'     # Gaia XP Giants
    }
    return mapping.get(priority, 'GIANTS')

def create_weave_compatible_catalog(pisces_file, output_file):
    """Create WEAVE-compatible catalog from Pisces targets."""
    
    print("Loading Pisces catalog...")
    with fits.open(pisces_file) as hdul:
        pisces_data = Table(hdul[1].data)
    
    print(f"Loaded {len(pisces_data)} Pisces targets")
    
    # Query Gaia data for all source IDs
    source_ids = pisces_data['SOURCE_ID'].data
    gaia_data = query_gaia_data(source_ids)
    
    if gaia_data is None:
        raise RuntimeError("Failed to retrieve Gaia data from wsdb")
    
    print(f"Retrieved Gaia data for {len(gaia_data)} targets")
    
    # Convert to astropy table and join
    gaia_table = Table.from_pandas(gaia_data)
    
    # Join on source_id - need to match column names  
    enriched = join(pisces_data, gaia_table, keys_left='SOURCE_ID', keys_right='source_id', join_type='left')
    
    print(f"Joined data: {len(enriched)} targets")
    
    # Create WEAVE-compatible columns
    n_targets = len(enriched)
    
    # Basic target information
    cname = np.array([f"PISCES_{i+1:06d}" for i in range(n_targets)], dtype='U20')
    targsrvy = np.full(n_targets, 'WS2025B2-028', dtype='U15')
    targcat = np.full(n_targets, 'WS2025B2-028.fits', dtype='U30')
    targid = np.array([f"PISCES_{enriched['SOURCE_ID'][i]}" for i in range(n_targets)], dtype='U30')
    targname = np.array([f"G{enriched['SOURCE_ID'][i]}" for i in range(n_targets)], dtype='U30')
    
    # Map priorities and target types
    targprog = np.array([map_priority_to_target_type(p) for p in enriched['PRIORITY']], dtype='U40')
    targprio = enriched['PRIORITY'].astype(np.float32)
    targuse = np.full(n_targets, 'T', dtype='U1')  # Target use flag
    targclass = np.copy(targprog).astype('U12')  # Same as TARGPROG for now
    
    # Observation parameters  
    progtemp = np.full(n_targets, '11331', dtype='U8')
    obstemp = np.full(n_targets, 'DACEB', dtype='U5')
    
    # Gaia information
    gaia_id = enriched['SOURCE_ID'].astype('U25')
    gaia_dr = np.full(n_targets, 'EDR3', dtype='U5')
    gaia_ra = enriched['ra'] if 'ra' in enriched.columns else enriched['RA']
    gaia_dec = enriched['dec'] if 'dec' in enriched.columns else enriched['DEC']
    gaia_epoch = np.full(n_targets, 2016.0, dtype=np.float32)
    
    # Handle missing values with defaults - use column access for astropy Table
    def safe_column(table, colname, default):
        if colname in table.columns:
            col = table[colname].data
            return np.where(np.isfinite(col), col, default)
        else:
            return np.full(len(table), default)
    
    gaia_pmra = safe_column(enriched, 'pmra', 0.0)
    gaia_pmra_err = np.full(n_targets, 999.0)  # Default error for missing data
    gaia_pmdec = safe_column(enriched, 'pmdec', 0.0)  
    gaia_pmdec_err = np.full(n_targets, 999.0)  # Default error for missing data
    gaia_paral = safe_column(enriched, 'parallax', -999.0)
    gaia_paral_err = np.full(n_targets, 999.0)  # Default error for missing data
    
    # Gaia magnitudes
    gaia_mag_g = safe_column(enriched, 'phot_g_mean_mag', 99.0)
    gaia_mag_g_err = np.full(n_targets, 9.0)  # Default error 
    gaia_mag_bp = safe_column(enriched, 'phot_bp_mean_mag', 99.0)
    gaia_mag_bp_err = np.full(n_targets, 9.0)  # Default error
    gaia_mag_rp = safe_column(enriched, 'phot_rp_mean_mag', 99.0)
    gaia_mag_rp_err = np.full(n_targets, 9.0)  # Default error
    
    # Galactic coordinates
    gaia_gal_lat = safe_column(enriched, 'gal_lat', -999.0)
    gaia_gal_long = safe_column(enriched, 'gal_long', -999.0)
    
    # HEALPix (compute from coordinates)
    import healpy as hp
    nside = 128
    healpix = hp.ang2pix(nside, np.radians(90 - gaia_dec), np.radians(gaia_ra))
    
    # Targeting information
    ga_targbits = np.zeros(n_targets, dtype=np.int64)  # Placeholder
    ga_targdate = np.full(n_targets, '2025-01-01', dtype='U10')
    ga_targrev = np.full(n_targets, 1, dtype=np.int32)
    
    # Extinction - use placeholder values since ebv column not available
    ebvcat = np.full(n_targets, 'GAIA', dtype='U10')
    ebvcat_dec = gaia_dec
    ebvcat_dr = np.full(n_targets, 'EDR3', dtype='U10')
    ebvcat_ebv = np.full(n_targets, -999.0)  # Placeholder
    ebvcat_ebv_err = np.full(n_targets, -999.0)  # Not available
    ebvcat_ra = gaia_ra
    
    # Create output table with all required columns
    output_table = Table()
    
    # Add all columns in order matching WEAVE template
    # First add the required ID columns that the tool expects
    output_table['SOURCE_ID'] = enriched['SOURCE_ID'].astype(np.int64)  # Keep original SOURCE_ID
    output_table['PS1_ID'] = enriched['PS1_ID'].astype(np.int64) if 'PS1_ID' in enriched.columns else np.zeros(n_targets, dtype=np.int64)
    output_table['RA'] = gaia_ra  # Add RA for compatibility
    output_table['DEC'] = gaia_dec  # Add DEC for compatibility 
    output_table['PRIORITY'] = targprio  # Keep PRIORITY column
    output_table['GAIA_REV_ID'] = np.full(n_targets, 3, dtype=np.int64)  # Gaia DR3 = 3
    
    output_table['CNAME'] = cname
    output_table['TARGSRVY'] = targsrvy
    output_table['TARGPROG'] = targprog
    output_table['TARGCAT'] = targcat
    output_table['TARGID'] = targid
    output_table['TARGNAME'] = targname
    output_table['TARGPRIO'] = targprio
    output_table['TARGUSE'] = targuse
    output_table['TARGCLASS'] = targclass
    output_table['PROGTEMP'] = progtemp
    output_table['OBSTEMP'] = obstemp
    output_table['GAIA_ID'] = gaia_id
    output_table['GAIA_DR'] = gaia_dr
    output_table['GAIA_RA'] = gaia_ra
    output_table['GAIA_DEC'] = gaia_dec
    output_table['GAIA_EPOCH'] = gaia_epoch
    output_table['GAIA_PMRA'] = gaia_pmra
    output_table['GAIA_PMRA_ERR'] = gaia_pmra_err
    output_table['GAIA_PMDEC'] = gaia_pmdec
    output_table['GAIA_PMDEC_ERR'] = gaia_pmdec_err
    output_table['GAIA_PARAL'] = gaia_paral
    output_table['GAIA_PARAL_ERR'] = gaia_paral_err
    output_table['GAIA_MAG_G'] = gaia_mag_g
    output_table['GAIA_MAG_G_ERR'] = gaia_mag_g_err
    output_table['GAIA_MAG_BP'] = gaia_mag_bp
    output_table['GAIA_MAG_BP_ERR'] = gaia_mag_bp_err
    output_table['GAIA_MAG_RP'] = gaia_mag_rp
    output_table['GAIA_MAG_RP_ERR'] = gaia_mag_rp_err
    output_table['GAIA_GAL_LAT'] = gaia_gal_lat
    output_table['GAIA_GAL_LONG'] = gaia_gal_long
    output_table['HEALPIX'] = healpix
    output_table['GA_TARGBITS'] = ga_targbits
    output_table['GA_TARGDATE'] = ga_targdate
    output_table['GA_TARGREV'] = ga_targrev
    output_table['EBVCAT'] = ebvcat
    output_table['EBVCAT_DEC'] = ebvcat_dec
    output_table['EBVCAT_DR'] = ebvcat_dr
    output_table['EBVCAT_EBV'] = ebvcat_ebv
    output_table['EBVCAT_EBV_ERR'] = ebvcat_ebv_err
    output_table['EBVCAT_RA'] = ebvcat_ra
    
    # Add placeholder columns for remaining WEAVE template requirements
    # (MAG_G, MAG_R, MAG_I, PS/SDSS cross-matches, etc.)
    placeholder_float = np.full(n_targets, -999.0)
    placeholder_str = np.full(n_targets, '', dtype='U10')
    
    output_table['MAG_G'] = placeholder_float.copy()
    output_table['MAG_G_ERR'] = placeholder_float.copy()
    output_table['MAG_R'] = placeholder_float.copy() 
    output_table['MAG_R_ERR'] = placeholder_float.copy()
    output_table['MAG_I'] = placeholder_float.copy()
    output_table['MAG_I_ERR'] = placeholder_float.copy()
    
    print(f"Created enriched catalog with {len(output_table)} targets and {len(output_table.columns)} columns")
    
    # Write to FITS file
    output_table.write(output_file, format='fits', overwrite=True)
    print(f"Saved enriched catalog to: {output_file}")
    
    return output_table

if __name__ == "__main__":
    pisces_file = "/Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/weave_pisces/Pisces_021023.fits"
    output_file = "/Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/weave_pisces/Pisces_enriched_021023.fits"
    
    enriched_catalog = create_weave_compatible_catalog(pisces_file, output_file)