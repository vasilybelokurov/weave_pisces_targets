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
weave_pisces_targets/
├── notebooks/
│   ├── pisces_targets_updated.ipynb      # Main target selection pipeline
│   ├── pisces_targets.ipynb              # Original implementation
│   ├── gaia_edr3_2mass_query.ipynb      # Expanded RGB catalog creation
│   ├── pisces_rrlyrae_map.ipynb          # RR Lyrae analysis for context
│   ├── cv_coord.py                       # Coordinate conversion utilities
│   ├── pisces_targets.py                 # Python target selection module
│   └── sphere_rotate.py                  # Spherical coordinate rotation
├── idl_reference/
│   ├── gaia_edr3_targets_pisces_proposal.pro  # Original IDL pipeline
│   └── low_parallax_query_2mass_ebv.pro       # Initial data query
├── docs/
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
├── Pisces_021023.fits                   # Final target catalog (5,191 targets)
├── gaiadr3_rrab_data.fits              # RR Lyrae reference catalog
├── CLAUDE.md                            # Development workflow
├── PROJECT_DESCRIPTION.md               # Project overview
├── LICENSE                              # MIT license
└── README.md                            # This file
```

## Target Categories

The pipeline selects four main categories of stellar tracers:

### 1. Clean RGB Giants
- **Selection**: `parallax < 0.2 mas, RPM < 18.0, (BP-RP)₀ > 1.0, G₀ < 18.0, dist > 50 kpc`
- **Count**: 531 targets in survey region
- **Priority**: 10 (highest)
- **Purpose**: High-confidence distant giants for kinematic analysis

### 2. Relaxed RGB Giants  
- **Selection**: `parallax < 0.25 mas, RPM < 18.5, (BP-RP)₀ > 1.0, G₀ < 18.5, dist > 30 kpc`
- **Count**: 2,467 additional targets
- **Priority**: 9
- **Purpose**: Expanded sample for increased statistics

### 3. Blue Horizontal Branch Stars
- **Selection**: `PS1 g > 18.0` (distance proxy)
- **Count**: 1,402 targets  
- **Priority**: 8
- **Purpose**: Alternative distance tracers and population indicators

### 4. Gaia XP Giants
- **Selection**: `ebv < 0.1, logg < 2.0, Teff < 5000, [M/H] < -1.5, G < 17.0, RPM < 17.0`
- **Count**: 791 targets
- **Priority**: 7
- **Purpose**: Spectroscopically confirmed giants with known stellar parameters

**Total**: 5,191 targets across 35 WEAVE pointings (2° diameter fields)

## Survey Parameters

- **Survey Region**: RA = [-45°, +15°], Dec = [-20°, 0°]
- **Area**: 1,176 deg² (spherical)
- **Instrument**: WEAVE multi-object spectrograph (WHT)
- **Field of View**: 2° diameter per pointing
- **Coverage**: 35 pointings (9.4% of survey area)
- **Expected**: ~14 targets per pointing

## Installation & Usage

### Requirements
- Python 3.8+
- Jupyter Notebook
- Standard scientific Python stack (numpy, pandas, matplotlib, astropy)
- Optional: sqlutilpy for database queries

### Quick Start

1. **Clone the repository**:
```bash
git clone https://github.com/vasilybelokurov/weave_pisces_targets.git
cd weave_pisces_targets
```

2. **Install dependencies**:
```bash
pip install -r requirements.txt  # Will be added
```

3. **Run main pipeline**:
```bash
jupyter notebook notebooks/pisces_targets_updated.ipynb
```

### Data Requirements

The pipeline uses several external catalogs:
- **Gaia EDR3**: RGB giant photometry and astrometry
- **2MASS**: Near-infrared photometry for color-magnitude selection  
- **PS1**: BHB star identification
- **Gaia XP**: Stellar parameter estimates from BP/RP spectra

Most data is accessed via web URLs or can be queried from astronomical databases.

## Key Features

### Scientific
- **Multi-catalog integration**: Gaia EDR3, 2MASS, PS1, XP spectra
- **Rigorous selection criteria**: Parallax, proper motion, color-magnitude cuts
- **Distance estimation**: Isochrone fitting using multiple stellar populations
- **Extinction correction**: Wavelength-dependent corrections for accurate photometry

### Technical  
- **Modular design**: Configuration-driven parameter system
- **Efficient implementation**: Vectorized operations and boolean masking
- **Comprehensive visualization**: Sky maps, selection diagnostics, parameter optimization
- **Quality control**: Finite value checking, statistical validation
- **Professional output**: FITS files with full metadata

### Analysis Tools
- **RR Lyrae context mapping**: Survey region validation using Pisces Plume tracers
- **Parameter space optimization**: 2D histograms comparing target/comparison regions
- **Selection diagnostics**: Color-magnitude, reduced proper motion analysis
- **Coverage analysis**: Expected targets per pointing based on surface density

## Results

The pipeline produces a target catalog (`Pisces_YYMMDD.fits`) with:
- **5,191 stellar targets** across 4 categories
- **Optimized for WEAVE**: Magnitude limits, spatial distribution
- **Science priorities**: Hierarchical ranking system
- **Complete metadata**: Gaia source IDs, coordinates, selection flags

Expected observational outcomes:
- **~490 targets per night** (5 nights total)
- **Magnitude range**: G = 13-18.5, suitable for WEAVE
- **Radial velocity precision**: ~1 km/s for kinematic analysis
- **Metallicity estimates**: [Fe/H] precision ~0.2 dex

## Documentation

Detailed documentation is provided in the `docs/` directory:

- **[Notebook Comparison](docs/PISCES_TARGETS_COMPARISON.md)**: Evolution from original to updated pipeline
- **[RGB Selection](docs/EXPANDED_RG_SELECTION.md)**: Methodology for expanded giant catalog
- **[Development Workflow](docs/CLAUDE.md)**: Coding standards and practices

## Citation

If you use this target selection pipeline in your research, please cite:

```bibtex
@misc{weave_pisces_targets,
  author = {Belokurov, Vasily and {WEAVE Pisces Team}},
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

- **Lead**: Dr. Vasily Belokurov (Institute of Astronomy, Cambridge)
- **Email**: vasily@ast.cam.ac.uk  
- **WEAVE Consortium**: [https://ingconfluence.ing.iac.es/confluence/display/WEAV](https://ingconfluence.ing.iac.es/confluence/display/WEAV)

## Acknowledgments

This work is supported by the WEAVE consortium and uses data from:
- European Space Agency (ESA) mission Gaia
- Two Micron All Sky Survey (2MASS)  
- Panoramic Survey Telescope and Rapid Response System (Pan-STARRS1)
- Isaac Newton Group of Telescopes, La Palma