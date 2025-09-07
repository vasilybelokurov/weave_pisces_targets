# WEAVE Pisces Target Selection: Notebook Comparison

## Overview

This document provides a detailed comparison between two Jupyter notebooks for WEAVE Pisces target selection:

- **`pisces_targets.ipynb`** (Original) - Initial implementation
- **`pisces_targets_updated.ipynb`** (Updated) - Refactored and enhanced version

## Structure and Architecture Comparison

### **pisces_targets.ipynb (Original)**

**Structure**: Linear, procedural approach with 37 cells
- Simple imports and data loading
- Sequential target selection with basic cuts
- Basic plotting and visualization
- Minimal modularity and reusability

**Key Characteristics**:
- **Ad-hoc approach**: Selection criteria hardcoded throughout cells
- **Minimal documentation**: Limited comments and explanations
- **Basic plotting**: Simple scatter plots, limited diagnostic capabilities
- **Limited flexibility**: Parameters scattered across cells, difficult to modify
- **Simple target packing**: Basic FITS file creation

### **pisces_targets_updated.ipynb (Updated)**

**Structure**: Modular, configuration-driven approach with 16 cells
- Centralized configuration system
- Modular functions and helper utilities  
- Comprehensive plotting and analysis framework
- Enhanced documentation and structure

**Key Characteristics**:
- **Configuration-driven**: Centralized parameter dictionaries
- **Modular design**: Reusable functions and utilities
- **Comprehensive analysis**: Advanced plotting, optimization tools, summary tables
- **Professional structure**: Clear separation of concerns, enhanced documentation
- **Flexible architecture**: Easy to modify parameters and extend functionality

## Detailed Cell-by-Cell Comparison

### **Data Loading and Setup**

| **Original** | **Updated** |
|--------------|-------------|
| **Cells 0-5**: Basic imports, manual URL loading | **Cell 0**: Professional imports with plotting style |
| Hardcoded URLs, simple data loading | **Cell 2**: Flexible loading with local/remote options |
| Manual RA conversion with explicit loops | **Cell 3**: Efficient helper functions |

**Key Improvements**:
- **Centralized imports** with consistent styling
- **Flexible data sources** (local expanded catalog vs remote)
- **Efficient coordinate transformations** using vectorized operations

### **Configuration and Parameters**

| **Original** | **Updated** |
|--------------|-------------|
| **Cells 13, 17**: Scattered parameter definitions | **Cell 1**: Centralized configuration dictionaries |
| ```python | ```python |
| rpm_cut=18. | rgb_params = { |
| par_cut=0.2 | 'clean': {'rpm_cut': 18.0, 'par_cut': 0.2, ...}, |
| col_cut=1. | 'relaxed': {'rpm_cut': 18.5, 'par_cut': 0.25, ...} |
| ``` | } |

**Key Improvements**:
- **Hierarchical parameter organization** by target type
- **Priority system** integrated into configuration
- **Easy parameter modification** without code changes
- **Multiple selection modes** (clean vs relaxed)

### **Target Selection Implementation**

| **Original** | **Updated** |
|--------------|-------------|
| **Cells 14-15**: Manual intersection operations | **Cell 5**: Vectorized boolean masks |
| ```python | ```python |
| wpar=np.where(dg['parallax']<par_cut) | def rgb_mask(kind='clean'): |
| wcol=np.where(dg['col']>col_cut) |     p = rgb_params[kind] |
| wtar1=reduce(np.intersect1d,(wpar,wcol,...)) |     return (dg['parallax'] < p['par_cut']) & ... |
| ``` | ``` |

**Key Improvements**:
- **Boolean mask approach** - more efficient and readable
- **Parameterized functions** for different selection criteria
- **Global vs field separation** handled systematically
- **Automatic additional relaxed** target computation

### **Analysis and Optimization**

| **Original** | **Updated** |
|--------------|-------------|
| **Missing**: No optimization tools | **Cell 4**: RR Lyrae density analysis with survey region comparison |
| **Missing**: No parameter space analysis | **Cell 4**: Multi-panel RRab density maps with old/new survey regions |
| Basic plotting only | **Cell 6**: RGB Selection Optimizer with 2D parameter space analysis |

**Key Improvements**:
- **RR Lyrae context maps** showing survey region in relation to Pisces Plume
- **Parameter space optimization** with difference plots (NEW - COMP regions)
- **Multi-dimensional analysis** across color-magnitude-RPM space
- **Threshold overlay visualization** showing selection boundaries

### **Plotting and Visualization**

| **Original** | **Updated** |
|--------------|-------------|
| **Cells 16-25**: Basic scatter plots | **Cell 8**: Professional full-sky density maps |
| Simple matplotlib calls | **Cells 9-11**: Comprehensive plotting suite |
| No density maps | Log-scaled 2D histograms with proper colorbars |
| Limited sky coverage visualization | Survey region overlays on all plots |

**Key Improvements**:
- **2D density maps** instead of scatter plots for large datasets
- **Logarithmic scaling** for better dynamic range
- **Consistent survey region overlays** across all plots
- **Galactic and Sagittarius coordinate** transformations
- **Professional styling** with proper axis labels and colorbars

### **Summary and Reporting**

| **Original** | **Updated** |
|--------------|-------------|
| **Cell 34**: Basic count printing | **Cell 6**: Comprehensive summary table |
| ```python | ```python |
| print('Total number of targets ',ntar1+ntar12+ntar2) | # Spherical area calculations |
| ``` | # Expected targets per pointing |
| **Cell 35**: Added later, basic table | # Efficiency calculations |
| | # Professional formatted output |

**Key Improvements**:
- **Spherical geometry calculations** for accurate areas
- **Expected targets per pointing** based on surface density
- **Efficiency analysis** (fraction of global sample in survey box)
- **Professional table formatting** with comprehensive statistics
- **Coverage analysis** relating to WEAVE instrument specifications

### **Data Output and Storage**

| **Original** | **Updated** |
|--------------|-------------|
| **Cells 29-33**: Manual target packing | **Cell 12**: Vectorized target assembly |
| Repetitive code for each target type | Clean function-based approach |
| Basic FITS creation | Enhanced FITS with metadata |

**Key Improvements**:
- **Modular packing functions** eliminating code duplication
- **Vectorized operations** for better performance
- **Comprehensive metadata** in output files
- **Consistent naming conventions** and file organization

## Target Selection Differences

### **Target Categories**

| **Category** | **Original** | **Updated** | **Change** |
|--------------|--------------|-------------|------------|
| **Clean RGB** | 531 targets | 531 targets | **Identical** |
| **Relaxed RGB** | 2,467 additional | 2,467 additional | **Identical** |
| **BHB** | 1,402 targets | 1,402 targets | **Identical** |
| **XP Giants** | Not included | 791 targets | **New category** |
| **Total** | 4,400 targets | 5,191 targets | **+791 (+18%)** |

### **Enhanced Target Categories**

The updated notebook introduces **Gaia XP giants** as a new target category:
- **Selection criteria**: `ebv<0.1, logg<2.0, Teff<5000, [M/H]<-1.5, G<17.0, RPM<17.0`
- **Priority**: 7 (lowest priority)
- **Contribution**: 791 additional targets (2.1 per pointing expected)

## Methodological Improvements

### **1. Parameter Space Optimization**
- **NEW vs COMP region analysis** for validating selection criteria
- **2D histogram comparisons** showing target density differences
- **Threshold visualization** overlaying selection cuts on parameter spaces

### **2. Survey Strategy Enhancement**
- **RR Lyrae context mapping** to validate survey region placement
- **Multiple survey region proposals** with comparative analysis
- **Spherical geometry** for accurate area and coverage calculations

### **3. Scientific Rigor**
- **Extinction correction** consistently applied across all catalogs
- **Coordinate system transformations** (Galactic, Sagittarius)
- **Distance estimation** using multiple methods
- **Quality control** with finite value checks and error handling

### **4. Computational Efficiency**
- **Boolean masks** instead of index-based selections
- **Vectorized operations** throughout
- **Memory-efficient** data loading with memmap
- **Reduced computational complexity** in target selection

## Code Quality and Maintainability

### **Original Notebook Issues**:
- **Hardcoded parameters** scattered across multiple cells
- **Repetitive code** for similar operations
- **Limited error handling** and validation
- **Inconsistent naming** conventions
- **Poor documentation** and comments

### **Updated Notebook Improvements**:
- **Configuration-driven** design with centralized parameters
- **DRY principle** applied with reusable functions
- **Comprehensive error handling** and input validation
- **Consistent coding standards** and naming conventions
- **Extensive documentation** with clear cell descriptions
- **Professional structure** suitable for production use

## Performance and Scalability

### **Original**:
- **Linear complexity** with manual intersections
- **Memory inefficient** with intermediate arrays
- **Limited to predefined** target types

### **Updated**:
- **Optimized boolean operations** for better performance
- **Memory-mapped file access** for large catalogs
- **Modular design** easily extensible to new target types
- **Scalable architecture** for larger datasets

## Conclusion

The **`pisces_targets_updated.ipynb`** represents a significant advancement over the original notebook:

1. **Scientific Enhancement**: More rigorous methodology with comprehensive analysis tools
2. **Engineering Excellence**: Professional code structure with improved maintainability  
3. **Operational Improvement**: Better target yield (+18%) with new XP giant category
4. **User Experience**: Clear documentation, modular design, and easy parameter modification
5. **Research Value**: Enhanced visualization and optimization tools for survey planning

The updated notebook transforms a functional but basic target selection script into a comprehensive, professional-grade survey planning tool suitable for large astronomical projects like WEAVE Pisces.