# Data dictionary

| Path | Contents | Used by |
| --- | --- | --- |
| `data/geno_with_pheno.csv` | Yeast cross genotypes and morphology phenotypes. | Figures 4 and 5B |
| `data/Correlations.Rfile` | Precomputed genetic-correlation table for the yeast morphology analysis. | Figures 4 and 5B |
| `data/genotype_yeast.rdata` | Yeast genotype matrix (`X`). | Figures 4 and 5B |
| `data/LD_pruning.rdata` | Marker names excluded by the LD-pruning analysis (`out_LD`). | Figure 4 |
| `data/merged_correlation.Rfile` | Genetic-correlation summaries across geldanamycin concentrations. | Figure 5B |
The yeast morphology data derive from Geiler-Samerotte et al. (2020). The
precomputed `.Rfile` and `.rdata` objects are analysis inputs, not raw-source
archives.
