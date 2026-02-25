# Code for Main Figures

Code and data to reproduce the main figures in:

**"Consistent and idiosyncratic pleiotropy: dissecting the genetic architecture of trait correlations through a genetic cross"**

Haoran Cai, Kerry Geiler-Samerotte, David L. Des Marais

## Repository Structure

```
public_code/
├── README.md
├── toolbox.r            # Shared utility functions (theme, simulation, bootstrap)
├── Fig3.r               # Figure 3: rG vs rD under simulated genetic architectures
├── Fig5.r               # Figure 5: Empirical rG vs rD, IQR outlier removal, Table 1
├── Fig6.r               # Figure 6: rD predicts rG stability under perturbation
└── data/
    ├── geno_with_pheno.csv         # Yeast cross genotype + phenotype data
    ├── Correlations.Rfile          # Pre-computed genetic correlations
    ├── genotype_yeast.rdata        # Yeast genotype matrix
    ├── LD_pruning.rdata            # LD-pruned locus list
    └── merged_correlation.Rfile    # Correlations across GdA concentrations
```

## Requirements

Install the following R packages before running the scripts:

```r
install.packages(c(
  "tidyverse", "qtl", "ggpubr", "viridis", "ggExtra", "corrr",
  "ggthemes", "confintr", "broom", "jtools", "bbplot", "patchwork",
  "PhenotypeSimulator", "MASS", "car", "readxl", "cowplot",
  "ggrepel", "diffcor", "reshape2", "gridExtra", "here",
  "preprocessCore", "Biobase", "broom.mixed"
))

# Bioconductor packages (if not already installed):
# BiocManager::install(c("preprocessCore", "Biobase"))
```

## Usage

All scripts use the [`here`](https://here.r-lib.org/) package for portable path management. To run any script:

1. Open R/RStudio
2. Set the working directory to this folder, or simply open the `.r` file — `here` will auto-detect the project root
3. Source the script:

```r
source("Fig3.r")    # generates Fig3_1.pdf, Fig3_2.pdf
source("Fig5.r")    # generates Fig5_1.pdf, Fig5_2.pdf, Horizontal_trait.csv
source("Fig6.r")    # generates Fig6_1.pdf, Fig6_2.png, Fig6_SI.pdf
```

## Data Sources

- **Yeast morphology data**: From [Geiler-Samerotte et al. (2020)](https://doi.org/10.7554/eLife.58564)
- **Brassica floral trait data**: From [Brock et al. (2010)](https://doi.org/10.1111/j.1558-5646.2010.01089.x)
