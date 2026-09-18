# Published-version analysis code

This repository contains the selected analyses supporting the published version
of _Consistent and idiosyncratic pleiotropy: dissecting the genetic architecture
of trait correlations through a genetic cross_ by Haoran Cai, Kerry
Geiler-Samerotte, and David L. Des Marais.

## Scope

The repository deliberately reproduces only the analyses retained for the
published article:

| Script | Publication content | Products |
| --- | --- | --- |
| `scripts/figure_4_trimmed_rd.R` | IQR-trimmed polygenic background-correlation analysis. | Figure 4 analytical panels and candidate-pair table |
| `scripts/figure_5b_perturbation_stability.R` | Figure 5B perturbation-stability analysis and supporting null analyses. | Figure 5B and supporting files |

The following are intentionally excluded: the retired legacy simulation,
published Figure 3, Figure 1's conceptual illustration, and the workflow
diagram in Figure 5A. Figure 5A is a manuscript graphic; this repository
reproduces its quantitative panel (5B), not the diagram.

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
  "car", "cowplot", "here", "tidygraph", "igraph", "broom.mixed"
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

Each script creates `results/` if necessary and writes its outputs there. The
empirical scripts use fixed random seeds where stochastic estimates are made.

## Data provenance

The yeast morphology inputs derive from
[Geiler-Samerotte et al. (2020)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000836).
See [the data dictionary](docs/data_dictionary.md) for each supplied object and
the scripts that use it.

## License

No license has been selected yet. Reuse permission is therefore not granted by
this repository. Add a license before creating the public release.
