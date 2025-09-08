# WEAVE Pisces Target Selection

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Jupyter](https://img.shields.io/badge/jupyter-notebook-orange.svg)](https://jupyter.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Target selection pipeline for the WEAVE Pisces survey, focusing on RGB giants and other stellar tracers in the outer Galactic halo.

## Overview

This repository contains the complete target selection pipeline for the WEAVE Pisces project, designed to study the outer halo structure and the Pisces Plume overdensity. The pipeline selects multiple categories of stellar tracers optimized for spectroscopic follow-up with the WEAVE multi-object spectrograph.

## Scientific Context

The Pisces Plume is a prominent stellar overdensity in the outer Galactic halo (40-110 kpc) first discovered through RR Lyrae stars. Recent work by [Belokurov et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488L..47B) suggests it may be part of the Magellanic wake - a gravitational wake induced by the Large Magellanic Cloud's infall into the Milky Way.

### Survey Goals
- **Map outer halo structure** using RGB giants and blue horizontal branch stars
- **Characterize the Pisces Plume** through multi-object spectroscopy
- **Study Magellanic Cloud influence** on Galactic halo dynamics
- **Measure kinematics and metallicities** of distant halo tracers

## Repository Structure

```
weave_pisces/
├── notebooks/
│   ├── pisces_targets_updated.ipynb      # Main target selection pipeline
│   ├── pisces_targets.ipynb              # Original implementation
│   ├── gaia_edr3_2mass_query.ipynb      # Expanded RGB catalog creation
│   ├── pisces_rrlyrae_map.ipynb          # RR Lyrae analysis for context
│   ├── cv_coord.py                       # Coordinate conversion utilities
│   ├── pisces_targets.py                 # Python target selection module
│   └── sphere_rotate.py                  # Spherical coordinate rotation
├── weave_tool/                           # Self-contained WEAVE catalog generation
│   ├── enrich_pisces_catalog.py          # Database enrichment script
│   ├── weave_pisces_script.sh            # Main generation script
│   ├── weave_fits_writeout/              # Adapted WEAVE chervin tool
│   ├── weave_pisces_env/                 # Local Python environment
│   ├── requirements.txt                  # Dependencies
│   └── docs/                             # Tool documentation
├── idl_reference/
│   ├── gaia_edr3_targets_pisces_proposal.pro  # Original IDL pipeline
│   └── low_parallax_query_2mass_ebv.pro       # Initial data query
├── docs/
│   ├── LATEST_CHANGES.md                # Major updates and 3.25x target increase analysis
│   ├── PROJECT_DESCRIPTION.md           # Complete scientific overview
│   ├── RG_OPTIMIZER.md                  # RGB parameter optimization tool guide
│   ├── PISCES_TARGETS_COMPARISON.md     # Detailed notebook comparison
│   ├── EXPANDED_RG_SELECTION.md         # RGB selection methodology
│   └── CLAUDE.md                        # Development workflow
├── proposals/
│   ├── 2025B-2026A_CCI_ITP_Pisces.pdf  # WEAVE observing proposal
│   ├── 2025B-2026A_CCI_ITP_Pisces.docx # Original proposal document
│   ├── 2025B-2026A_CCI_ITP_Pisces.xml  # Proposal metadata
│   ├── Belokurov_etal_2019.pdf          # Pisces Plume discovery paper
│   └── Trager_Scott_WEAVEstatus-20250525-short.pdf  # WEAVE status
├── plots/
│   ├── pisces_all_targets_combined.png  # Combined target distribution
│   ├── pisces_*_fullsky.png             # Full-sky target maps
│   ├── pisces_*_zoom.png                # Survey region zoom plots
│   ├── pisces_clean_rgb_sgr_*.png       # Sagittarius coordinate projections
│   ├── rg_selection_optimizer*.png      # Parameter optimization plots
│   └── pisces_targets_distribution.png  # Target category breakdown
├── Pisces_021023.fits                   # Target catalog (14,952 targets) - Version 1
├── Pisces_080925.fits                   # Target catalog (14,952 targets) - Latest version
├── Pisces_WEAVE_enhanced_021023.fits    # Enhanced WEAVE catalog - Version 1  
├── Pisces_WEAVE_enhanced_080925.fits    # Enhanced WEAVE catalog - Latest version
├── gaiadr3_rrab_data.fits              # RR Lyrae reference catalog
├── LICENSE                              # MIT license
└── README.md                            # This file
```

## Target Categories

The pipeline selects five main categories of stellar tracers:

### 1. Clean RGB Giants
- **Selection**: `parallax < 0.2 mas, RPM < 19.5, (BP-RP)₀ > 1.0, G₀ < 19.5, dist > 25 kpc`
- **Count**: 5,471 targets in optimized survey region
- **Priority**: 10 (highest)
- **Purpose**: High-confidence distant giants for kinematic analysis

### 2. Additional Relaxed RGB Giants  
- **Selection**: `parallax < 0.25 mas, RPM < 20.0, (BP-RP)₀ > 0.9, G₀ < 20.0, dist > 15 kpc`
- **Count**: 8,028 additional targets
- **Priority**: 9
- **Purpose**: Expanded sample with relaxed cuts for comprehensive halo coverage

### 3. Blue Horizontal Branch Stars
- **Selection**: `PS1 g > 18.0` (distance proxy)
- **Count**: 794 targets  
- **Priority**: 8
- **Purpose**: Alternative distance tracers and population indicators

### 4. Gaia XP Giants
- **Selection**: `ebv < 0.1, logg < 2.0, Teff < 5000, [M/H] < -1.5, G < 17.0, RPM < 17.0`
- **Count**: 248 targets
- **Priority**: 7
- **Purpose**: Spectroscopically confirmed giants with known stellar parameters

### 5. RR Lyrae Variables (RRab)
- **Selection**: `distance_kpc > 30.0` (standard candles)
- **Count**: 411 targets
- **Priority**: 8
- **Purpose**: Independent distance tracers and kinematic validation

**Total**: 14,952 targets across 35 WEAVE pointings (2° diameter fields)

## Survey Parameters

- **Optimized Survey Region**: RA = [-20°, +10°], Dec = [-22°, +7°] (aligned with RR Lyrae overdensity)
- **Area**: 853 deg² (spherical)  
- **Instrument**: WEAVE multi-object spectrograph (WHT)
- **Field of View**: 2° diameter per pointing
- **Coverage**: 35 pointings (12.9% of survey area)
- **Expected**: ~55 targets per pointing

## WEAVE Catalog Generation Tool

This repository includes a **complete WEAVE-compatible catalog generation tool** that transforms basic Pisces target catalogs into full 78-column WEAVE fiber catalogs ready for spectroscopic observations.

### What the Tool Does

The WEAVE tool **enriches** your basic 6-column Pisces catalog:

**Input** (`Pisces_080925.fits` - latest version):
- Basic columns: `SOURCE_ID`, `RA`, `DEC`, `GAIA_REV_ID`, `PS1_ID`, `PRIORITY`
- 14,952 Pisces targets selected by priority

**Output** (`Pisces_WEAVE_enhanced_080925.fits`):
- **78 columns** with complete stellar data and WEAVE metadata
- **Cross-matched photometry** from Gaia EDR3, PS1, SDSS databases
- **Extinction corrections** using SFD dust maps
- **WEAVE survey format** ready for fiber allocation planning

### Key Enhancements Added

1. **Astrometry & Photometry**: Gaia EDR3 positions, proper motions, parallaxes, G/BP/RP magnitudes + errors
2. **Multi-band Photometry**: PS1 gri, SDSS ugriz magnitudes where available  
3. **WEAVE Metadata**: Survey identifiers, target classifications, program templates
4. **Quality Control**: Coordinate precision, extinction corrections, finite value validation
5. **Database Integration**: Live queries to wsdb.ast.cam.ac.uk for latest stellar data

### Usage

**One-command catalog generation:**
```bash
cd weave_tool
./weave_pisces_script.sh ../Pisces_080925.fits
# Processing ~15-30 minutes for 14,952 targets
# Output: WS2025B2-028-PISCES.fits → moved to project root as Pisces_WEAVE_enhanced_080925.fits
```

**Or from Jupyter notebook:**
```python
# Add this cell to pisces_targets_updated.ipynb (already integrated)
result = subprocess.run(
    ['./weave_pisces_script.sh', '../Pisces_080925.fits'],
    cwd='../weave_tool'
)
```

### Tool Components

All WEAVE-related files are organized in the **`weave_tool/`** directory:

- **`weave_pisces_script.sh`** - Main catalog generation script
- **`weave_fits_writeout/`** - Adapted WEAVE chervin tool for Pisces survey  
- **`weave_pisces_env/`** - Isolated Python environment with dependencies
- **`enrich_pisces_catalog.py`** - Database enrichment utilities
- **`requirements.txt`** - Complete dependency list with SFD dust maps
- **`docs/`** - Detailed setup and troubleshooting documentation

### Requirements

- **Database access**: wsdb.ast.cam.ac.uk (uses system PGUSER/PGPASSWORD)
- **Python 3.8+** with scientific computing stack
- **SFD dust maps**: Downloaded automatically on first run (~128 MB)
- **Processing time**: ~2 seconds per target (network dependent)

See **[WEAVE Tool Setup](weave_tool/docs/WEAVE_TOOL_SETUP.md)** for complete documentation.

## Installation & Usage

### Requirements

**For Target Selection:**
- Python 3.8+ with Jupyter Notebook
- Standard scientific Python stack (numpy, pandas, matplotlib, astropy)

**For WEAVE Catalog Generation:**
- Access to wsdb.ast.cam.ac.uk database (with PGUSER/PGPASSWORD configured)
- Internet connection for SFD dust map downloads (~128 MB, one-time)
- Processing time: ~15-30 minutes for full catalog

### Quick Start - Target Selection

1. **Clone the repository**:
```bash
git clone https://github.com/vasilybelokurov/weave_pisces_targets.git
cd weave_pisces
```

2. **Run target selection pipeline**:
```bash
jupyter notebook notebooks/pisces_targets_updated.ipynb
```

**Output**: `Pisces_080925.fits` (14,952 targets, 6 columns)

### Quick Start - WEAVE Catalog Generation

**After completing target selection above:**

1. **Generate enhanced WEAVE catalog**:
```bash
cd weave_tool
./weave_pisces_script.sh ../Pisces_080925.fits
```

**Output**: `Pisces_WEAVE_enhanced_080925.fits` (14,952 targets, 78 columns)

### Complete Workflow

```bash
# 1. Target Selection (Jupyter notebook)
jupyter notebook notebooks/pisces_targets_updated.ipynb
# → Creates Pisces_080925.fits

# 2. WEAVE Enhancement (automated)
cd weave_tool  
./weave_pisces_script.sh ../Pisces_080925.fits
# → Creates Pisces_WEAVE_enhanced_080925.fits

# 3. Ready for WEAVE observation planning!
```

### Advanced Usage

**Custom target catalogs**: Replace `Pisces_080925.fits` with your own 6-column catalog:
- Required columns: `SOURCE_ID`, `RA`, `DEC`, `GAIA_REV_ID`, `PS1_ID`, `PRIORITY`
- Format: FITS binary table

**Database configuration**: Tool uses system environment variables:
- `PGUSER` - Your wsdb username  
- `PGPASSWORD` - Your wsdb password
- `PGHOST` - Set to `wsdb.ast.cam.ac.uk` (automatic in script)

**Troubleshooting**: See `weave_tool/docs/` for detailed setup guides

### Data Requirements

The pipeline uses several external catalogs:
- **Gaia EDR3**: RGB giant photometry and astrometry
- **2MASS**: Near-infrared photometry for color-magnitude selection  
- **PS1**: BHB star identification
- **Gaia XP**: Stellar parameter estimates from BP/RP spectra

Most data is accessed via web URLs or can be queried from astronomical databases.

## Key Features

### Scientific
- **Multi-catalog integration**: Gaia EDR3, 2MASS, PS1, XP spectra, expanded RGB catalog
- **Optimized survey region**: Data-driven selection aligned with RR Lyrae overdensity
- **Five target populations**: RGB giants, BHB stars, XP giants, RR Lyrae variables
- **Expanded sample size**: 5.7M global RGB candidates (5.7x larger than standard catalog)
- **Distance estimation**: Isochrone fitting using multiple stellar populations
- **Extinction correction**: Wavelength-dependent corrections for accurate photometry

### Technical  
- **Modular design**: Configuration-driven parameter system
- **Efficient implementation**: Vectorized operations and boolean masking
- **Comprehensive visualization**: Sky maps, selection diagnostics, parameter optimization
- **Quality control**: Finite value checking, statistical validation
- **Professional output**: FITS files with full metadata

### Analysis Tools
- **RG Selection Optimizer**: Critical tool for data-driven parameter optimization using 2D histograms
- **RR Lyrae context mapping**: Survey region validation using Pisces Plume tracers  
- **Background subtraction**: Comparison with offset regions to reveal genuine overdensity
- **Selection diagnostics**: Color-magnitude, reduced proper motion analysis
- **Coverage analysis**: Expected targets per pointing based on surface density

## Results

The complete pipeline produces two complementary catalogs:

### Target Selection Catalog (`Pisces_080925.fits` - Latest)
- **14,952 stellar targets** across 5 categories in optimized survey region
- **6 essential columns**: SOURCE_ID, RA, DEC, GAIA_REV_ID, PS1_ID, PRIORITY
- **Science priorities**: Hierarchical ranking system (priorities 7-10)
- **Size**: ~710 KB, ready for target selection analysis

### Enhanced WEAVE Catalog (`Pisces_WEAVE_enhanced_080925.fits` - Latest)  
- **Same 14,952 targets** with complete observational metadata
- **78 detailed columns**: Full Gaia EDR3 astrometry, multi-band photometry, WEAVE survey parameters
- **Database-enriched**: Live cross-matches with Gaia, PS1, SDSS catalogs via wsdb.ast.cam.ac.uk
- **Extinction-corrected**: SFD dust map corrections for precise stellar parameters
- **WEAVE-ready**: Complete format compliance for fiber allocation and spectroscopic pipeline
- **Size**: ~10.5 MB, production-ready for WEAVE observations

### Previous Versions
- `Pisces_021023.fits` and `Pisces_WEAVE_enhanced_021023.fits` - Earlier versions maintained for reference

### Key Enhancements in WEAVE Catalog
- **Precision astrometry**: Gaia proper motions (±0.1 mas/yr), parallaxes (±0.1 mas)
- **Multi-band photometry**: G, BP, RP (Gaia) + gri (PS1) + ugriz (SDSS where available)
- **Survey metadata**: Target classifications, program templates (11331), observation modes
- **Quality assurance**: Error propagation, finite value validation, coordinate precision

### Expected Observational Outcomes
- **~1,495 targets per night** (5 nights total at 35 pointings)
- **Magnitude range**: G = 13-20, optimized for WEAVE spectroscopy
- **Target density**: 55.0 targets per pointing (approaches WEAVE fiber limits)
- **Sample size**: 3.25x larger than previous selection for robust statistics
- **Radial velocity precision**: ~1 km/s for kinematic analysis  
- **Metallicity estimates**: [Fe/H] precision ~0.2 dex
- **Survey efficiency**: Ready for immediate WEAVE observation planning

## Documentation

Detailed documentation is provided in the `docs/` directory:

- **[Latest Changes](docs/LATEST_CHANGES.md)**: Major updates and 3.25x target increase analysis
- **[Project Description](docs/PROJECT_DESCRIPTION.md)**: Complete scientific overview and technical specifications
- **[RG Optimizer](docs/RG_OPTIMIZER.md)**: Critical tool for RGB parameter optimization and validation
- **[Notebook Comparison](docs/PISCES_TARGETS_COMPARISON.md)**: Evolution from original to updated pipeline
- **[RGB Selection](docs/EXPANDED_RG_SELECTION.md)**: Methodology for expanded giant catalog
- **[Development Workflow](docs/CLAUDE.md)**: Coding standards and practices

## Citation

If you use this target selection pipeline in your research, please cite:

```bibtex
@misc{weave_pisces_targets,
  author = {Belokurov, Vasily and Deason, Alis J. and {WEAVE Pisces Team}},
  title = {WEAVE Pisces Target Selection Pipeline},
  year = {2025},
  url = {https://github.com/vasilybelokurov/weave_pisces_targets},
  note = {Target selection for the WEAVE Pisces outer halo survey}
}
```

## References

Key papers for scientific context:
- [Belokurov et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488L..47B) - "The Pisces Plume and the Magellanic wake"
- [Gaia Collaboration 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...649A...1G) - Gaia EDR3 data release
- [WEAVE Collaboration 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.5323D) - WEAVE survey overview

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Co-Principal Investigators:**
- **Dr. Vasily Belokurov** (Institute of Astronomy, Cambridge)
  - Email: vasily@ast.cam.ac.uk
- **Dr. Alis J. Deason** (Institute for Computational Cosmology, Durham University)  
  - Email: alis.j.deason@durham.ac.uk

**WEAVE Project**: [https://weave-project.atlassian.net/wiki/spaces/WEAVE/overview](https://weave-project.atlassian.net/wiki/spaces/WEAVE/overview)

## Acknowledgments

This work is supported by the WEAVE consortium and uses data from:
- European Space Agency (ESA) mission Gaia
- Two Micron All Sky Survey (2MASS)  
- Panoramic Survey Telescope and Rapid Response System (Pan-STARRS1)
- Isaac Newton Group of Telescopes, La Palma
