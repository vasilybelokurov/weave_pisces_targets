# RGB Giant Selection Optimizer

## Overview

The RGB giant selection optimizer is a critical diagnostic tool implemented in `pisces_targets_updated.ipynb` that enables data-driven optimization of target selection parameters. This tool compares the stellar density distributions between the target survey region and comparison regions to maximize signal-to-noise for detecting stellar substructure.

## Scientific Motivation

### The Problem
Traditional RGB giant selection relies on fixed cuts in color-magnitude and proper motion space. However, optimal parameters depend on:
- **Contamination levels** from nearby disk/halo populations
- **Signal strength** of the target overdensity
- **Background subtraction** effectiveness
- **Survey-specific magnitude limits** and observational constraints

### The Solution
The optimizer uses 2D histograms to visualize stellar density ratios between:
- **Target region**: Pisces survey box (RA=[-45°,+15°], Dec=[-20°,0°])
- **Comparison regions**: Offset boxes designed to sample background populations

## Implementation

### Current Parameters Being Optimized
1. **Parallax cut**: `parallax < 0.25 mas` (distance proxy)
2. **RPM limit**: `RPM < 18.5` (luminosity class separator)
3. **Color cut**: `(BP-RP)₀ > 1.0` (giant vs. dwarf separation)
4. **Magnitude limit**: `G₀ < 18.5` (observational constraint)
5. **Distance cut**: `distance > 30 kpc` (foreground rejection)

### Visualization Output
The optimizer generates density ratio maps showing:
- **Parameter space** (e.g., RPM vs. G magnitude)
- **Color-coded enhancement** (target/comparison density ratio)
- **Selection boundaries** overlaid on parameter distributions
- **Quantitative assessment** of target purity vs. completeness

## Key Findings

### Current Selection Efficiency
- **Global sample**: 721,246 RGB candidates (relaxed selection)
- **Survey region**: 2,998 targets (0.4% efficiency)
- **Background subtraction**: Reveals genuine Pisces overdensity

### Parameter Sensitivity
The optimizer reveals that:
1. **RPM cuts** are most effective for contamination rejection
2. **Distance estimates** show large systematic uncertainties
3. **Color cuts** provide robust giant/dwarf separation
4. **Magnitude limits** define the depth-completeness trade-off

## Critical Need for Exploration

### Why We Need to Revisit Selection
The current RGB selection was derived from the original IDL pipeline but may not be optimal for several reasons:

1. **Expanded catalog**: We now have access to the expanded RG catalog with improved photometry
2. **Better distance estimates**: Updated isochrone fitting and extinction corrections
3. **Refined background models**: More sophisticated comparison region sampling
4. **WEAVE constraints**: Instrument-specific magnitude and density requirements

### Recommended Parameter Space Exploration

#### High Priority
- **RPM vs. G magnitude**: Optimize the luminosity-distance relation
- **Color-magnitude**: Refine the RGB giant locus in extinction-corrected space
- **Distance cuts**: Balance foreground rejection vs. sample completeness

#### Medium Priority  
- **Parallax thresholds**: Test sensitivity to astrometric quality cuts
- **Extinction corrections**: Validate wavelength-dependent coefficients
- **Proper motion**: Explore kinematic selection for specific substructures

#### Low Priority
- **Multi-dimensional cuts**: Consider ML-based selection algorithms
- **Adaptive parameters**: Region-dependent selection optimization

## Usage Instructions

### Running the Optimizer
1. **Execute cells 1-6** in `pisces_targets_updated.ipynb` to load catalogs
2. **Modify parameters** in the configuration section (cell 3)
3. **Run the optimizer** (typically around cell 7-8, look for 2D histogram plots)
4. **Analyze output**: Examine density ratio maps and selection boundaries
5. **Iterate parameters** based on visual assessment and target counts

### Interpreting Results
- **Red/yellow regions**: Overdense in target vs. comparison (good signal)
- **Blue regions**: Underdense in target (potential contamination)
- **Contour lines**: Current selection boundaries
- **Target counts**: Balance between purity and sample size

## Integration with Target Selection

### Current Workflow
1. **Generate comparison regions** offset from target survey box
2. **Apply selection cuts** to both target and comparison samples  
3. **Compute density ratios** in 2D parameter space
4. **Visualize enhancement patterns** to validate selection
5. **Adjust parameters** based on optimization results

### Future Improvements
- **Automated parameter scanning**: Grid search over parameter combinations
- **Quantitative metrics**: Signal-to-noise ratio, completeness estimates
- **Statistical validation**: Bootstrap resampling, cross-validation
- **Integration with WEAVE**: Fiber allocation and observability constraints

## Files and Dependencies

### Primary Implementation
- **`pisces_targets_updated.ipynb`**: Main optimizer implementation
- **Helper modules**: `cv_coord.py`, `sphere_rotate.py` for coordinate transformations
- **Data catalogs**: RGB catalog (local or web-based), comparison samples

### Output Products
- **`plots/rg_selection_optimizer*.png`**: Visualization of parameter optimization
- **Selection diagnostics**: Target counts, efficiency metrics
- **Parameter logs**: Record of tested configurations

## Conclusion

The RGB optimizer is essential for validating and improving our target selection methodology. **Immediate exploration is needed** to ensure our final target list maximizes scientific return from the WEAVE Pisces survey. The tool provides both qualitative (visual) and quantitative (target counts) feedback to guide parameter optimization in a data-driven manner.

**Next steps**: Systematically explore the recommended parameter combinations, document the optimization process, and update selection criteria based on empirical validation using the optimizer output.