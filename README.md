# Consistent and idiosyncratic pleiotropy in shaping genetic correlations

This repository contains the selected analyses supporting the published version
of **Consistent and idiosyncratic pleiotropy in shaping genetic correlations**
by Haoran Cai, Kerry Geiler-Samerotte, and David L. Des Marais.

## Scope

The repository deliberately reproduces only the analyses retained for the
published article:

| Script | Publication content | Products |
| --- | --- | --- |
| `scripts/figure_4_trimmed_rd.R` | IQR-trimmed polygenic background-correlation analysis. | Figure 4 analytical panels and candidate-pair table |
| `scripts/figure_5b_perturbation_stability.R` | Figure 5B perturbation-stability analysis and supporting null analyses. | Figure 5B and supporting files |

## Layout

```text
R/                         shared analysis functions
scripts/                   figure-specific entry points
data/                      supplied empirical analysis inputs
results/                   regenerated files (ignored by Git)
docs/data_dictionary.md    input descriptions
```

## Setup

Use R 4.4 or newer. The Figure 4 and 5B analyses use the following packages:

```r
install.packages(c(
  "tidyverse", "qtl", "ggpubr", "viridis", "ggExtra", "corrr",
  "ggthemes", "confintr", "broom", "jtools", "bbplot", "patchwork",
  "car", "cowplot", "here", "tidygraph", "igraph", "broom.mixed",
  "LaplacesDemon", "MultiRNG", "VGAM"
))

install.packages("BiocManager")
BiocManager::install(c("preprocessCore", "Biobase"))
```

## Reproduce retained analyses

From the repository root, run the scripts in this order:

```r
source("scripts/figure_4_trimmed_rd.R")
source("scripts/figure_5b_perturbation_stability.R")
```
