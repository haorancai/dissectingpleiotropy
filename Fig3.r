# =============================================================================
# Figure 3: Simulation of rG vs rD under different genetic architectures
# =============================================================================
#
# This script generates Figure 3, which demonstrates the relationship between
# genetic correlation (rG) and polygenic background correlation (rD) under:
#   - Purely polygenic architecture (no HP)
#   - Architecture with large-effect horizontal pleiotropic (HP) loci
#   - Two kurtosis levels (gamma = 0.5 and gamma = 1)
#
# Data requirements:
#   - data/geno_with_pheno.csv  (yeast genotype + phenotype cross data)
#   - data/Correlations.Rfile   (pre-computed genetic correlations)
#
# Output:
#   - Fig3_1.pdf  (4-panel: polygenic vs HP, low vs high kurtosis)
#   - Fig3_2.pdf  (effect size distribution comparison)
# =============================================================================

# --- Load packages -----------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(bbplot)
library(PhenotypeSimulator)
library(qtl)
library(MASS)
library(patchwork)

# --- Set up paths (uses `here` for portability) ------------------------------
# install.packages("here")  # uncomment if not installed
library(here)
here::i_am("Fig3.r")
source(here("toolbox.r"))


# =============================================================================
# Panel: Effect size distribution under different kurtosis
# =============================================================================

NrSNP <- 226
mu <- c(0, 0)

n <- 20
set.seed(1010)
patch0 <- tibble(gamma = c(rep(0.5, n), rep(1, n))) %>%
    rowwise() %>%
    mutate(rau = runif(1, -1, 1)) %>%
    mutate(value = map2(gamma, rau, ~ draw.multivariate.laplace(NrSNP, 2, .x, mu, matrix(c(1, .y, .y, 1),
        ncol = 2
    ))[, 1])) %>%
    ungroup() %>%
    mutate(num = row_number()) %>%
    unnest(a = value) %>%
    mutate(gamma = factor(gamma)) %>%
    dplyr::select(gamma, a) %>%
    group_by(gamma) %>%
    mutate(a = scale(a)) %>%
    ggplot(aes(a, group = fct_rev(gamma), fill = fct_rev(gamma))) +
    geom_histogram(binwidth = .3, alpha = .5, position = "identity") +
    theme_Publication() +
    theme(
        legend.position = "top",
        axis.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")
    ) +
    xlab("Effect size") +
    ylab("Count") +
    scale_fill_manual(name = "Kurtosis", values = c("#798E87", "#9C964A"), labels = c("Low", "High")) +
    theme(aspect.ratio = 1 / 1)

# Color definitions
high <- "#9C964A"
low <- "#798E87"


# =============================================================================
# Panel: rG vs rD from empirical yeast data
# =============================================================================

NrSNP <- 226
sample_size <- 500

genotypes <- simulateGenotypes(
    N = sample_size, NrSNP = NrSNP,
    frequencies = 0.5,
    verbose = FALSE
)
SNPinfo <- genotypes$genotypes

# Load yeast cross data
cross <- read.cross(format = "csv", file = here("data", "geno_with_pheno.csv"), genotypes = c("A", "B"))
class(cross)[1] <- "dh"

additive_effect <- get_additive_effect(cross)
load(here("data", "Correlations.Rfile"))

df <- Correlations %>%
    as_tibble() %>%
    dplyr::select(Phen1, Phen2) %>%
    mutate(a = map2(Phen1, Phen2, ~ tibble(
        a1 = additive_effect %>% select(.x) %>% pull(),
        a2 = additive_effect %>% select(.y) %>% pull()
    ))) %>%
    mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
    unnest(cor_a) %>%
    rowwise() %>%
    mutate(cor_g = cor(SNPinfo %*% cbind(a$a1, a$a2))[2, 1])

df %>%
    ggplot(aes(cor_a, cor_g)) +
    geom_point(alpha = .1, color = "Black") +
    xlab("rD") +
    ylab("rG") +
    theme_Publication() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank()
    ) +
    theme(aspect.ratio = 1 / 1) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    geom_vline(xintercept = 0, size = 1, colour = "#333333")

ggsave("Fig3_s1.pdf")


# =============================================================================
# Simulations: Polygenic (no HP) architecture
# =============================================================================

NrSNP <- 226
sample_size <- 500

genotypes <- simulateGenotypes(
    N = sample_size, NrSNP = NrSNP,
    frequencies = 0.5,
    verbose = FALSE
)
SNPinfo <- genotypes$genotypes

# High kurtosis (gamma = 1), polygenic
df <- c()
for (i in c(1:1000)) {
    df <- rbind(df, simulation_for_a_given_population_exponential(
        NrSNP = NrSNP, SNPinfo = SNPinfo, rau = NULL,
        effect_size = FALSE, mode = "polygenic", gamma = 1
    ))
}

p1 <- df %>%
    as_tibble() %>%
    rename(genetic_cor = V1, developmental_cor = V2) %>%
    ggplot(aes(developmental_cor, genetic_cor)) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    geom_vline(xintercept = 0, size = 1, colour = "#333333") +
    geom_point(alpha = .5, color = low) +
    xlab("rD") +
    ylab("rG") +
    theme_Publication() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank()
    ) +
    theme(aspect.ratio = 1 / 1)

# Low kurtosis (gamma = 0.5), polygenic
df <- c()
for (i in c(1:1000)) {
    df <- rbind(df, simulation_for_a_given_population_exponential(
        NrSNP = NrSNP, SNPinfo = SNPinfo, rau = NULL,
        effect_size = FALSE, mode = "polygenic", gamma = 0.5
    ))
}

p2 <- df %>%
    as_tibble() %>%
    rename(genetic_cor = V1, developmental_cor = V2) %>%
    ggplot(aes(developmental_cor, genetic_cor)) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    geom_vline(xintercept = 0, size = 1, colour = "#333333") +
    geom_point(alpha = .65, color = high) +
    xlab("rD") +
    ylab("rG") +
    theme_Publication() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank()
    ) +
    theme(aspect.ratio = 1 / 1)


# =============================================================================
# Simulations: Major-effect (HP) architecture
# =============================================================================

# High kurtosis (gamma = 1), major-effect HP loci
df <- c()
for (i in c(1:1000)) {
    df <- rbind(df, simulation_for_a_given_population_exponential(
        NrSNP = NrSNP, SNPinfo = SNPinfo, rau = NULL,
        effect_size = FALSE, mode = "major", gamma = 1
    ))
}

p3 <- df %>%
    as_tibble() %>%
    rename(genetic_cor = V1, developmental_cor = V2) %>%
    ggplot(aes(developmental_cor, genetic_cor)) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    geom_vline(xintercept = 0, size = 1, colour = "#333333") +
    geom_point(alpha = .5, color = low) +
    xlab("rD") +
    ylab("rG") +
    theme_Publication() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank()
    ) +
    theme(aspect.ratio = 1 / 1)

# Low kurtosis (gamma = 0.5), major-effect HP loci
df <- c()
for (i in c(1:1000)) {
    df <- rbind(df, simulation_for_a_given_population_exponential(
        NrSNP = NrSNP, SNPinfo = SNPinfo, rau = NULL,
        effect_size = FALSE, mode = "major", gamma = 0.5
    ))
}

p4 <- df %>%
    as_tibble() %>%
    rename(genetic_cor = V1, developmental_cor = V2) %>%
    ggplot(aes(developmental_cor, genetic_cor)) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    geom_vline(xintercept = 0, size = 1, colour = "#333333") +
    geom_point(alpha = .65, color = high) +
    xlab("rD") +
    ylab("rG") +
    theme_Publication() +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank()
    ) +
    theme(aspect.ratio = 1 / 1)


# =============================================================================
# Assemble and save
# =============================================================================

patch1 <- p1 + p2 + grid::textGrob("No HP") + plot_layout(widths = c(2, 2, 1))
patch2 <- p3 + p4 + grid::textGrob("HP") + plot_layout(widths = c(2, 2, 1))

Fig3_1 <- patch1 / patch2
ggsave(Fig3_1, filename = "Fig3_1.pdf")
ggsave(patch0, filename = "Fig3_2.pdf", width = 3, height = 3)
