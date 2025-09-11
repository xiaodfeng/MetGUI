# MetGUI: An R Package for Comprehensive LC-MS Metabolomics Data Analysis

Detailed materials can be accessed via 10.6084/m9.figshare.30090904.

MetGUI is an **open-source R package** with a graphical user interface (GUI) designed for end-to-end pre-processing and analysis of LC-MS-based metabolomics data. Developed for researchers with limited programming experience, MetGUI integrates state-of-the-art algorithms into a user-friendly workflow that handles raw data conversion, peak detection, metabolite identification, and statistical analysis.

## Key Features
- **Local Installation**: Process sensitive data without cloud uploads
- **Comprehensive Workflow**: Six integrated modules covering all analysis steps
- **Advanced Algorithms**:
  - Dynamic binning for accurate peak detection
  - Optimized scaling factors for metabolite identification
  - m/z-dependent peak matching with dynamic tolerance
- **Visualization Tools**: Interactive diagnostic plots at each processing stage
- **FAIR Compliance**: Parameter logging for reproducible research

## Installation
1. Ensure **R (≥ 4.2.3)** is installed
2. Install required dependencies:
```r
install.packages(c("shiny", "MsnBase", "xcms", "CAMERA", "MetaboAnalystR", "msPurity"))
```
3. Install MetGUI from GitHub:
```r
devtools::install_github("xiaodfeng/MetGUI")
```

## Quick Start
```r
library(MetGUI)
launch_metgui()  # Starts the interactive interface
```

## Workflow Modules
### 1. Data Transformation
- Convert vendor formats → open standards (mzML)
- Visualize TIC/BPC chromatograms
- Project setup and parameter management

### 2. Peak Detection & Quantification
- Dynamic binning algorithm accounts for:
  - Mass analyzer type (Orbitrap, Q-TOF, etc.)
  - Instrument-specific resolving power
  - m/z-dependent peak location uncertainty

### 3. Spectra Pre-processing
- Retention time alignment (OBI-warp algorithm)
- Peak grouping across samples
- Missing value imputation
- Diagnostic EIC visualizations

### 4. Annotation
- Isotope pattern recognition and filtering
- Adduct annotation and neutral mass calculation
- Duplicate peak removal
- Interactive scatter plots (RT vs. m/z)

### 5. Metabolite Identification
- Spectral library matching with:
  - Dynamic tolerance peak matching
  - Optimized cosine similarity scoring
- FDR-controlled identification
- Target-decoy visualizations
- Supports experimental and in-silico libraries

### 6. Statistics & Pathway Analysis
- Data normalization (6 methods)
- Volcano plots and multivariate analysis (PCA, PLS-DA)
- Pathway enrichment via Mummichog algorithm
- Automated report generation

## Case Study: Hypoxic vs. Normoxic K562 Cells
MetGUI was validated using an LC-MS/MS dataset of K562 cells under normoxic (21% O₂) and hypoxic (1% O₂) conditions:



*Key findings:*
- Identified 453 metabolites at 5% FDR
- Detected arginine upregulation in hypoxia
- Pathway analysis revealed altered glutathione metabolism
- Reduced TCA cycle activity under low oxygen

## Documentation
- vignettes/MetGUI manual.Rmd
- vignettes/DYNAMIC_BINNING.md

## System Requirements
- OS: Windows 10/11, Linux, macOS
- RAM: 16GB+ recommended for large datasets
- Storage: SSD recommended for raw file handling

## License
This project is licensed under the **GNU General Public License v3.0** - see LICENSE for details.

## Citation
If you use MetGUI in your research, please cite:  
Feng X, Cunningham A, et al. *MetGUI: an R package for metabolomics LC-MS data analysis*. (Manuscript in preparation)

## Contact
- **Lead Developer**: Dr. Xiaodong Feng (fxd@jclab.ac.cn)
- **GitHub Issues**: https://github.com/xiaodfeng/MetGUI/issues
