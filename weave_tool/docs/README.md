# WEAVE Tool for Pisces Survey

## Quick Start

Generate WEAVE-compatible catalog from Pisces targets:

```bash
cd /path/to/weave_pisces
./weave_pisces_script.sh Pisces_021023.fits
```

**Output**: `WS2025B2-028-PISCES.fits` (ready for WEAVE observation planning)

## Files

- **`WEAVE_TOOL_SETUP.md`** - Complete setup and usage documentation
- **`../enrich_pisces_catalog.py`** - Database enrichment script
- **`../weave_pisces_script.sh`** - Main processing script  
- **`../requirements.txt`** - Python dependencies
- **`../weave_pisces_env/`** - Local environment

## What This Tool Does

1. **Enriches** basic Pisces catalog with Gaia EDR3 data from wsdb
2. **Converts** to WEAVE-compatible format (78 columns)
3. **Maps** Pisces priorities (7-10) to WEAVE target types (GIANTS/BHB/RRLYRAE)
4. **Generates** final catalog ready for WEAVE fiber allocation

## Target Mapping

| Pisces Priority | Targets | WEAVE Type | Description |
|-----------------|---------|------------|-------------|
| 10 | 5,471 | GIANTS | Clean RGB Giants |
| 9  | 8,028 | GIANTS | Relaxed RGB Giants |
| 8  | 1,205 | BHB/RRLYRAE | BHB + RR Lyrae |
| 7  | 248   | GIANTS | Gaia XP Giants |

**Total**: 14,952 Pisces targets → WEAVE-compatible catalog