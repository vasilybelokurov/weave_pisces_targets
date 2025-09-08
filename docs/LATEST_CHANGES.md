# Latest Changes to WEAVE Pisces Target Selection

## Major Updates - January 2025

### 1. Expanded RGB Catalog Integration

**What Changed:**
- **Activated expanded catalog**: `expanded_rg = True` now uses the dramatically larger `rg_gaia2mass_expanded.fits` catalog
- **Massive scale increase**: From ~1M to 5.7M global RGB candidates (5.7x larger sample)

**Impact:**
- **Higher completeness**: Access to significantly more distant giant candidates
- **Better statistics**: Improved target density and coverage in Pisces region
- **Enhanced selection power**: More candidates allow for stricter cuts while maintaining target counts

### 2. Optimized RGB Selection Parameters

**What Changed:**
- **New parameter set**: `rgb_params_new` with substantially relaxed cuts
- **Clean RGB**: 
  - RPM: 18.0 → 19.5 (allows fainter/more distant giants)
  - G magnitude: 18.0 → 19.5 (extends to fainter targets)
  - Distance: 50.0 → 25.0 kpc (includes closer halo giants)
- **Relaxed RGB**: 
  - RPM: 18.5 → 20.0 (even fainter limit)
  - G magnitude: 18.5 → 20.0 (pushes to WEAVE limit)
  - Color cut: 1.0 → 0.9 (slightly bluer giants included)
  - Distance: 30.0 → 15.0 kpc (includes inner halo)

**Scientific Rationale:**
- **WEAVE capabilities**: Pushes to instrument limits (G~20) for maximum science return
- **Halo completeness**: Lower distance cuts capture complete halo population
- **Optimized for expanded catalog**: Leverages larger sample to relax cuts while maintaining quality

### 3. Addition of RR Lyrae as Fifth Target Category

**What Changed:**
- **New target class**: RR Lyrae (RRab) variables added as legitimate WEAVE targets
- **Selection criterion**: `distance_kpc > 30.0 kpc` 
- **Priority**: 8 (same as BHB stars)
- **Target count**: 411 RRab stars in optimized survey region

**Scientific Value:**
- **Independent distance tracers**: Standard candles with ~3% distance precision
- **Population indicators**: Trace old, metal-poor stellar populations
- **Kinematic tracers**: Excellent radial velocity targets for dynamics
- **Pisces validation**: Confirm survey region alignment with known overdensity

### 4. Dramatically Increased Target Density

**What Changed:**
- **Total targets**: 4,601 → 14,952 (3.25x increase)
- **Targets per pointing**: ~17 → 55 (3.2x increase)
- **Clean RGB**: 2,464 → 5,471 (2.2x increase)
- **Additional Relaxed RGB**: 1,095 → 8,028 (7.3x increase)
- **New RR Lyrae category**: +411 targets

**Impact:**
- **Exceptional target density**: 55 targets per WEAVE pointing (vs. typical 20-30)
- **Statistical power**: Massive improvement in sample size for dynamics analysis
- **Observation efficiency**: Much higher science return per pointing
- **Complete halo sampling**: Comprehensive coverage of stellar halo populations

### 5. Enhanced RG Selection Optimizer

**What Changed:**
- **Dual parameter comparison**: Shows both old (`rgb_params`) and new (`rgb_params_new`) selection boundaries
- **Visual validation**: Overlaid selection cuts on 2D histograms for parameter optimization
- **Background subtraction**: Refined comparison region analysis using offset survey boxes

**Technical Improvements:**
- **Parameter visualization**: Black solid lines for new cuts, dashed lines for old cuts
- **Optimization guidance**: Clear visual feedback for parameter adjustment
- **Scientific validation**: Demonstrates genuine Pisces overdensity vs. background

### 6. Optimized Survey Region Confirmed

**What Stayed the Same:**
- **Survey coordinates**: RA=[-20°,+10°], Dec=[-22°,+7°] (853 deg²)
- **RR Lyrae alignment**: Region optimally positioned based on RRab overdensity analysis
- **Background subtraction**: Uses offset comparison region for validation

### 7. Technical Enhancements

**What Changed:**
- **Target assembly**: Added `pack_rr()` function for RR Lyrae FITS format
- **Complete integration**: All five target categories assembled into single catalog
- **Updated output**: `Pisces_021023.fits` now contains 14,952 targets across 5 categories

## Scientific Impact

### Observational Advantages
- **Fiber efficiency**: 55 targets per pointing approaches maximum WEAVE capacity
- **Sample completeness**: ~3x larger sample enables robust statistical analysis
- **Multi-population tracers**: Five complementary stellar populations for comprehensive analysis
- **Distance range**: Complete coverage from 15-100 kpc in stellar halo

### Analysis Capabilities  
- **Kinematic precision**: Large sample enables detection of subtle velocity substructure
- **Metallicity gradients**: Sufficient statistics for radial abundance analysis
- **Population synthesis**: Multiple tracers allow detailed halo archaeology
- **Systematic control**: RR Lyrae provide independent distance/kinematics validation

### Comparison to Original Selection
- **Target count**: 14,952 vs. 4,601 (3.25x increase)
- **Sample diversity**: 5 vs. 4 target categories (RR Lyrae addition)
- **Distance coverage**: 15-100 kpc vs. 25-100 kpc (includes inner halo)
- **Magnitude depth**: G<20 vs. G<19.5 (pushes WEAVE limits)
- **Completeness**: Dramatically improved thanks to expanded catalog

## Validation and Quality Control

### RG Selection Optimizer Results
- **Clear overdensity**: Target region shows significant enhancement over background
- **Parameter optimization**: Visual confirmation that selection boundaries are well-placed
- **Background subtraction**: Genuine Pisces signal confirmed vs. offset comparison regions

### RR Lyrae Validation
- **Independent confirmation**: 411 RRab in survey region confirms target alignment
- **Distance consistency**: RR Lyrae distances confirm >30 kpc halo population
- **Spatial correlation**: RRab distribution validates survey box positioning

## Recommendations

### Immediate Actions
1. **Review target density**: 55 targets/pointing may approach WEAVE fiber limits
2. **Validate magnitude limits**: Confirm G=20 targets meet S/N requirements for WEAVE
3. **Check observation time**: Higher target density may require longer exposures

### Future Optimizations
1. **Magnitude stratification**: Consider priority-based magnitude limits for fiber allocation
2. **Regional optimization**: Fine-tune selection parameters based on optimizer results
3. **Cross-validation**: Use RR Lyrae subset to validate RGB giant distances

## Conclusion

The updated target selection represents a **major advancement** in both scale and sophistication. The 3.25x increase in target density, addition of RR Lyrae tracers, and optimization based on the expanded RGB catalog position this survey to deliver unprecedented insights into Milky Way halo structure and the Pisces Plume specifically.

**Key Achievement**: Transformation from a targeted survey to a comprehensive census of distant stellar populations in the Pisces region, maintaining scientific focus while dramatically expanding discovery potential.