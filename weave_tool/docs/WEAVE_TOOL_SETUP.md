# WEAVE Tool Setup for Pisces Survey

## Summary

Successfully adapted the existing `weave_chervin` tool to process Pisces target catalogs and generate WEAVE-compatible fiber catalogs.

## What Was Done

### 1. Environment Setup
- **Local virtual environment**: `weave_pisces_env/` 
- **Requirements file**: `requirements.txt` with all dependencies
- **Clean installation**: No conflicts with system packages

### 2. Tool Fixes Applied
- **NumPy 2.0 compatibility**: Changed `np.string_` → `np.bytes_` in `merge_cat.py`
- **Database connection**: Uses `WSDB_HOST=wsdb.ast.cam.ac.uk` with existing user credentials
- **Dust maps**: Downloaded SFD dust maps to local environment
- **Template creation**: Added `WS2025B2-028-PISCES_CatalogueTemplate.fits`

### 3. Pisces Catalog Enrichment
- **Created**: `enrich_pisces_catalog.py` - Enriches basic Pisces catalog with wsdb data
- **Input**: `Pisces_021023.fits` (14,952 targets, basic columns)
- **Output**: `Pisces_enriched_021023.fits` (15,184 targets, 52 columns)
- **Added required columns**: `SOURCE_ID`, `PS1_ID`, `RA`, `DEC`, `GAIA_REV_ID`, etc.

### 4. Pisces-Specific Configuration
- **Target mapping**: Created `pisces_target.yaml` with Pisces priority→WEAVE object type mapping
- **General config**: `pisces_general.yaml` with survey-specific parameters
- **Script**: `weave_pisces_script.sh` - One-command catalog generation

## Files Created/Modified

### New Files
```
weave_pisces_env/                    # Local virtual environment
requirements.txt                     # Python dependencies  
enrich_pisces_catalog.py             # Catalog enrichment script
Pisces_enriched_021023.fits          # Enriched target catalog (15,184 targets)
weave_pisces_script.sh               # Main processing script
WEAVE_TOOL_SETUP.md                  # This documentation

# Configuration files
weave_fits_writeout/py/weave_chervin/conf/pisces_target.yaml
weave_fits_writeout/py/weave_chervin/conf/pisces_general.yaml
weave_fits_writeout/py/weave_chervin/data/WS2025B2-028-PISCES_CatalogueTemplate.fits
```

### Modified Files
```
weave_fits_writeout/py/weave_chervin/merge_cat.py    # NumPy 2.0 fix
```

## Usage

### Quick Start
```bash
# Process Pisces targets to WEAVE format
cd /path/to/weave_pisces
./weave_pisces_script.sh Pisces_enriched_021023.fits
```

### Step-by-Step
```bash
# 1. Activate environment
source weave_pisces_env/bin/activate

# 2. Enrich basic catalog (if needed)
python enrich_pisces_catalog.py

# 3. Generate WEAVE catalog  
./weave_pisces_script.sh Pisces_enriched_021023.fits

# Output: WS2025B2-028-PISCES.fits (WEAVE-compatible)
```

## Target Mapping

| Pisces Priority | Target Count | WEAVE Object Type | Description |
|-----------------|--------------|-------------------|-------------|
| 10 | 5,471 | GIANTS | Clean RGB Giants |
| 9  | 8,028 | GIANTS | Additional Relaxed RGB Giants |
| 8  | 1,205 | BHB/RRLYRAE | BHB Stars + RR Lyrae |
| 7  | 248   | GIANTS | Gaia XP Giants |

**Total**: 14,952 → 15,184 targets (after enrichment and database joins)

## Database Requirements

- **Host**: `wsdb.ast.cam.ac.uk` (via `WSDB_HOST` environment variable)
- **Credentials**: Uses existing system `PGUSER`/`PGPASSWORD`
- **Required access**: Gaia EDR3, PS1, SDSS cross-match tables

## Output Format

- **File**: `WS2025B2-028-PISCES.fits`
- **Structure**: 2 HDUs (primary + binary table)
- **Columns**: 78 columns matching WEAVE template
- **Survey ID**: `WS2025B2-028-PISCES`
- **Program**: `11331` (Low-resolution spectroscopy)

## Validation

### Test Results
- **Test catalog**: 100 targets processed successfully
- **Database queries**: Working (PS1, Gaia, SDSS crossmatches)
- **Dust correction**: SFD maps integrated
- **Output format**: Validates against WEAVE template

### Key Metrics
- **Processing time**: ~2-3 minutes for 100 targets
- **Success rate**: 100% (all targets processed)
- **Output size**: ~1.2 MB per 100 targets

## Next Steps

1. **Full catalog processing**: Run on complete 15,184 target catalog
2. **Quality validation**: Check output against WEAVE requirements
3. **Fiber allocation**: Submit to WEAVE observation planning
4. **Documentation**: Update main project README with WEAVE workflow

## Dependencies

See `requirements.txt` for complete list. Key packages:
- `astropy>=5.0` - FITS file handling
- `sqlutilpy>=0.25` - Database queries  
- `healpy>=1.16` - HEALPix operations
- `dustmaps>=1.0` - Extinction corrections
- `weave_chervin` - WEAVE catalog generation (local install)

## Notes

- Tool preserves original `weave_chervin` functionality
- Local environment isolates dependencies
- All database credentials use existing system settings
- Compatible with WEAVE observation planning system