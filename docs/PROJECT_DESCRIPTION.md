# WEAVE Pisces Project Description

## Executive Summary
**Project:** Comprehensive spectroscopic survey of the Pisces Overdensity using the WEAVE spectrograph to disentangle three major galactic accretion events in the Milky Way's outer halo (50-110 kpc).

**Key Science Goals:**
- First spectroscopic detection of the Magellanic wake through radial velocity signatures
- Characterize the oldest tidal debris from the Sagittarius dwarf galaxy
- Map the elusive Magellanic stellar stream  
- Constrain the dynamical history and masses of both the Magellanic Clouds and Sagittarius dwarf

**Observational Parameters:** 35 WEAVE pointings × 60 min exposures = 5 nights total (2 dark + 3 grey)
**Target density:** 5,191 stellar targets across 4 categories with G = 13-18.5 mag, distances 30-110 kpc

## Scientific Context

### The Pisces Overdensity: A Stellar Crossroads
The Pisces Overdensity is a stellar enhancement in the outer Galactic halo spanning RA ≈ [-45°, +15°] and Dec ≈ [-20°, 0°]. Discovered through excess RR Lyrae and blue horizontal branch stars, this structure extends from 40-110 kpc with a steep distance gradient aligned with the Magellanic Stream.

Recent theoretical models suggest this region represents a unique convergence of three major stellar substructures:

1. **Magellanic Wake:** Dark matter and stellar debris focused by the LMC's gravitational passage
2. **Sagittarius Trailing Arm:** Oldest tidal debris from the Sgr dwarf's disruption  
3. **Magellanic Stellar Stream:** Long-sought stellar counterpart to the HI Magellanic Stream

### Critical Observational Need
Current spectroscopic data are extremely limited—only a handful of radial velocities exist for this region. The area lies outside the WEAVE Galactic Archaeology survey footprint (Dec < 0°), making this dedicated program essential for understanding these structures.

## Technical Specifications

### WEAVE Instrument Capabilities
- **Telescope:** William Herschel Telescope (4.2m)
- **Field of view:** 2° diameter
- **Fibers:** 960 (plate A) / 940 (plate B) 
- **Fiber size:** 1.3"
- **Resolution:** Low-resolution mode R ≈ 5000 (blue: 4150, red: 5950)
- **Wavelength coverage:** 3600-9490 Å
- **Current status:** Science verification ongoing (mid-2025), MOS mode available 2025B

### Target Selection Strategy

**Primary Targets (Priority 10):** Clean distant red giants (531 targets)
- Parallax < 0.2 mas, RPM < 18.0, (BP-RP)₀ > 1.0, G₀ < 18.0, distance > 50 kpc

**Secondary Targets (Priority 9):** Additional relaxed RGB giants (2,467 targets)  
- Parallax < 0.25 mas, RPM < 18.5, (BP-RP)₀ > 1.0, G₀ < 18.5, distance > 30 kpc

**Tertiary Targets (Priority 8):** Blue horizontal branch stars (1,402 targets)
- PS1 g > 18.0 mag for distance confirmation

**Quaternary Targets (Priority 7):** Gaia XP giants (791 targets)
- ebv < 0.1, logg < 2.0, Teff < 5000K, [M/H] < -1.5, G < 17.0, RPM < 17.0

**Total sample:** 5,191 targets providing ~14 targets per WEAVE pointing

### Expected Precision
- **Radial velocities:** ~1 km/s precision for kinematic analysis (G = 13-18.5)
- **Metallicities:** [Fe/H] precision ~0.2 dex  
- **α-element abundances:** Standard WEAVE pipeline output

## Science Impact & Legacy Value

### Immediate Science Returns
- **First detection of dynamical friction signature** in stellar halo through wake velocity structure
- **Sgr dwarf initial mass constraints** from oldest trailing debris metallicity/kinematics  
- **Magellanic orbital reconstruction** from stellar stream properties
- **Comprehensive chemo-dynamical atlas** of outer halo substructure

### Broader Implications
- Quantify recent massive accretion impact on Milky Way equilibrium
- Test theoretical models of satellite-halo interactions  
- Provide critical constraints for Galaxy formation simulations
- Enable accurate correction for Magellanic perturbations in future surveys

### Legacy Dataset
This program will deliver the most comprehensive spectroscopic mapping of distant halo stars (>30 kpc) to date, with multi-catalog target selection spanning RGB giants, BHB stars, and spectroscopically-confirmed giants. The resulting dataset provides a valuable resource for stellar halo studies and serves as ground truth for theoretical models of galactic accretion.

### Repository and Data Products
The complete target selection pipeline, including 4 Jupyter notebooks, parameter optimization tools, and visualization outputs, is publicly available at: https://github.com/vasilybelokurov/weave_pisces_targets

**Key outputs:**
- `Pisces_021023.fits`: Final target catalog (5,191 targets) with priorities and metadata
- `gaiadr3_rrab_data.fits`: RR Lyrae reference catalog for survey validation  
- Comprehensive documentation and methodology comparison with original IDL pipeline

## Team Expertise
**Principal Investigator:** Prof. Alis Deason (Durham University)  
**Key Collaborators:** V. Belokurov (Cambridge), G. Battaglia (IAC), S. Koposov (Edinburgh), D. Erkal (Surrey)

The team includes WEAVE Galactic Archaeology survey leads, target selection specialists, and world experts in Milky Way stellar halo dynamics, ensuring optimal survey execution and data analysis.