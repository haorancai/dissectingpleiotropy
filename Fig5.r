# =============================================================================
# Figure 5: Empirical rG vs rD with IQR outlier removal
# =============================================================================
#
# This script generates Figure 5 from the manuscript.
#   - Scatter plots of rG vs effect size correlation (default, LD-pruned,
#     outlier-corrected), plus density and deviation plots.
#
# Data requirements:
#   - data/geno_with_pheno.csv
#   - data/Correlations.Rfile
#   - data/genotype_yeast.rdata
#   - data/LD_pruning.rdata
#
# Output:
#   - Fig5_1.pdf, Fig5_2.pdf (Figure 5 panels)
#   - Horizontal_trait.csv (Table 1: candidate HP trait pairs)
# =============================================================================

# --- Load packages -----------------------------------------------------------
library(tidyverse)
library(qtl)
library(ggpubr)
library(viridis)
library(preprocessCore)
library(Biobase)
library(ggExtra)
library(corrr)
library(ggthemes)
library(confintr)
library(broom)
library(jtools)
library(bbplot)
library(patchwork)

# --- Set up paths ------------------------------------------------------------
library(here)
here::i_am("Fig5.r")
source(here("toolbox.r"))

# load data -------------------------------------------------------------------


cross <- read.cross(format = "csv", file = here("data", "geno_with_pheno.csv"), genotypes = c("A", "B"))
class(cross)[1] <- "dh"
load(here("data", "LD_pruning.rdata"))


# get_additive_effect -----------------------------------------------------

set.seed(1010)
additive_effect <- get_additive_effect(cross)


# Default ----------------------------------------------------------------

load(here("data", "Correlations.Rfile"))

# LOD <- LOD[!(locus_name %in% out_LD),]
Correlations <-
  Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull()
  ))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a)

Default <- Correlations$cor_a
Default_df <- Correlations %>%
  as_tibble() %>%
  # mutate(r_2_a =  ifelse(cor_a > 0, cor_a ** 2, -cor_a ** 2)) %>%
  mutate(Between = as.numeric(Between))

p3 <- Default_df %>%
  # mutate(number_of_large_qtl = a) %>%
  # filter(abs(Within) > 0.3) %>%
  #  mutate(Between = abs(Between)) %>%
  # filter(Between > 0.5) %>%
  ggplot(aes(cor_a, Between)) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  # scale_color_viridis(name = 'Environmental \n correlation',option="magma") +
  geom_point(alpha = .4, color = "#809F8C") +
  # theme_minimal()+#geom_smooth( method = lm, se = T)+
  # stat_function(fun = square,
  #              n = 200, color = 'steelblue4')+
  stat_function(
    fun = linear,
    n = 200, color = "black", size = 1.5
  ) +
  # scale_colour_gradient2()+
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2')+
  # stat_cor(method = "pearson")
  xlab("Effect Size Correlation \n (All variants)") +
  ylab("rG") +
  # ylim(0,1)+
  # gghighlight::gghighlight(label_key = Level, label_params = list(size = 8))+#, MINUTES %in% 60:90)+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  ) +
  # stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication() +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  geom_vline(xintercept = 0, size = 1, colour = "#333333") +
  # theme(title = element_text(face = 'bold'),
  #       #plot.title = element_text(color="black", size=5, face="bold"),
  #       axis.title.x = element_text(color="black", size=22),
  #       axis.title.y = element_text(color="black", size=22),
  #       legend.text = element_text(color = 'black', size = 10),
  #       legend.position = 'top',
  #       legend.title = element_text(color = 'black', size = 10))+
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  ) +
  theme(aspect.ratio = 1 / 1)


# LD ----------------------------------------------------------------

load(here("data", "Correlations.Rfile"))
load(here("data", "genotype_yeast.rdata"))
locus_name <- colnames(X)

# LOD <- LOD[!(locus_name %in% out_LD),]
Correlations <-
  Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  # sample_n(5) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    # lod1 = LOD %>% select(.x) %>% pull(),
    # lod2 = LOD %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    # mutate(order1 = row_number(-lod1),
    #        order2 = row_number(-lod2)) %>%
    # filter(order1 > 20 & order2 > 20) %>%
    filter(!(num %in% out_LD)))) %>%
  # select(-lod1,-lod2,-num))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a)

LD <- Correlations$cor_a

LD_df <- Correlations %>%
  as_tibble() %>%
  # mutate(r_2_a =  ifelse(cor_a > 0, cor_a ** 2, -cor_a ** 2)) %>%
  mutate(Between = as.numeric(Between))

p4 <- LD_df %>%
  ggplot(aes(cor_a, Between)) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  # scale_color_viridis(name = 'Environmental \n correlation',option="magma") +
  geom_point(alpha = .4, color = "#637095") +
  # theme_minimal()+#geom_smooth( method = lm, se = T)+
  # stat_function(fun = square,
  #              n = 200, color = 'steelblue4')+
  stat_function(
    fun = linear,
    n = 200, color = "black", size = 1.5
  ) +
  # scale_colour_gradient2()+
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2')+
  # stat_cor(method = "pearson")
  xlab("Effect Size Correlation \n (LD pruned variants)") +
  ylab("rG") +
  # ylim(0,1)+
  # gghighlight::gghighlight(label_key = Level, label_params = list(size = 8))+#, MINUTES %in% 60:90)+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  ) +
  # stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication() +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  geom_vline(xintercept = 0, size = 1, colour = "#333333") +
  # theme(title = element_text(face = 'bold'),
  #       #plot.title = element_text(color="black", size=5, face="bold"),
  #       axis.title.x = element_text(color="black", size=22),
  #       axis.title.y = element_text(color="black", size=22),
  #       legend.text = element_text(color = 'black', size = 10),
  #       legend.position = 'top',
  #       legend.title = element_text(color = 'black', size = 10))+
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  ) +
  theme(aspect.ratio = 1 / 1)


# Outlier ----------------------------------------------------------------


locus_name <- colnames(X)

load(here("data", "genotype_yeast.rdata"))
load(here("data", "Correlations.Rfile"))


# LOD <- LOD[!(locus_name %in% out_LD),]

Correlations <- Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  # sample_n(5) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == FALSE))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a)

outlier_corrected <- Correlations$cor_a


# bootstraping ------------------------------------------------------------


# D15.2_C: Bud nucleus brightness C
# D192_C
data_for_bootstrap <- Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  filter(Phen1 == "D15.2_C", Phen2 == "D192_C") %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  unnest(a)


set.seed(1010)
# [0.276, 0.547]
c(1:1000) %>%
  as_tibble() %>%
  mutate(a = map(value, ~ data_for_bootstrap %>%
    sample_n(226, replace = TRUE) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == FALSE))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a) %>%
  # summarize(sd =sd(cor_a))
  reframe(quants = quantile(cor_a, probs = c(0.025, 0.5, 0.975)))


# D117_C: D: nucleus center to cell center C
# D135_C: D: nucleus foci to cell center C
data_for_bootstrap <- Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  filter(Phen1 == "D117_C", Phen2 == "D135_C") %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  unnest(a)


set.seed(1010)
# [0.754, 0.888]
c(1:1000) %>%
  as_tibble() %>%
  mutate(a = map(value, ~ data_for_bootstrap %>%
    sample_n(226, replace = TRUE) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == FALSE))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a) %>%
  # summarize(sd =sd(cor_a))
  reframe(quants = quantile(cor_a, probs = c(0.025, 0.5, 0.975)))


# D128_C:D: bud nucleus foci to cell tip C
# D141_C:D: nucleus foci to cell hip C
data_for_bootstrap <- Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  filter(Phen1 == "D128_C", Phen2 == "D141_C") %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  unnest(a)


set.seed(1010)
# [0.981, 0.990]
c(1:1000) %>%
  as_tibble() %>%
  mutate(a = map(value, ~ data_for_bootstrap %>%
    sample_n(226, replace = TRUE) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == FALSE))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a) %>%
  # summarize(sd =sd(cor_a))
  reframe(quants = quantile(cor_a, probs = c(0.025, 0.5, 0.975)))


# D173_A:maxD: nucleus center to nucleus outline A
# D179_A:Minimum radius of nucleus A
data_for_bootstrap <- Correlations %>%
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within, Between) %>%
  filter(Phen1 == "D173_A", Phen2 == "D179_A") %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  unnest(a)


set.seed(1010)
# [0.963, 0.978]
c(1:1000) %>%
  as_tibble() %>%
  mutate(a = map(value, ~ data_for_bootstrap %>%
    sample_n(226, replace = TRUE) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == FALSE))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  unnest(cor_a) %>%
  # summarize(sd =sd(cor_a))
  reframe(quants = quantile(cor_a, probs = c(0.025, 0.5, 0.975)))


# a <- cbind(outlier_corrected_1.25,outlier_corrected_1.5) %>%
#   as_tibble() %>%
#   ggplot(aes(outlier_corrected_1.5, outlier_corrected_1.25))+geom_point()+theme_Publication()+
#   xlab('Outlier Corrected (cutoff = 1.5IQR)')+
#   ylab('Outlier Corrected (cutoff = 1.25IQR)')+
#   stat_cor(method = "pearson", label.x =-0.8, label.y = 1, r.accuracy = 0.0001)
#
#
#
# b <- cbind(outlier_corrected_1.75,outlier_corrected_1.5) %>%
#   as_tibble() %>%
#   ggplot(aes(outlier_corrected_1.5, outlier_corrected_1.75))+geom_point()+theme_Publication()+
#   xlab('Outlier Corrected (cutoff = 1.5IQR)')+
#   ylab('Outlier Corrected (cutoff = 1.75IQR)')+
#   stat_cor(method = "pearson", label.x =-0.8, label.y = 1, r.accuracy = 0.0001)
#
#
# library(patchwork)
# a + b
# ggsave(filename = 'Fig_S9.pdf')

# # theme_Publication()
# plot(outlier_corrected_1.25,outlier_corrected_1.5)
#
# plot(outlier_corrected_1.75,outlier_corrected_1.5)
#
# plot(outlier_corrected_1.25,outlier_corrected_1.75)
#
#
# cor(outlier_corrected_1.25,outlier_corrected_1.5)
#
# cor(outlier_corrected_1.75,outlier_corrected_1.5)
#
# cor(outlier_corrected_1.25,outlier_corrected_1.75)
#


outlier_corrected_df <- Correlations %>%
  as_tibble() %>%
  # mutate(r_2_a =  ifelse(cor_a > 0, cor_a ** 2, -cor_a ** 2)) %>%
  mutate(Between = as.numeric(Between))

p5 <- outlier_corrected_df %>%
  ggplot(aes(cor_a, Between)) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  # scale_color_viridis(name = 'Environmental \n correlation',option="magma") +
  geom_point(alpha = .4, color = "#AB6B6B") +
  # theme_minimal()+#geom_smooth( method = lm, se = T)+
  # stat_function(fun = square,
  #              n = 200, color = 'steelblue4')+
  stat_function(
    fun = linear,
    n = 200, color = "black", size = 1.5
  ) +
  # scale_colour_gradient2()+
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2')+
  # stat_cor(method = "pearson")
  xlab("Effect Size Correlation \n (Outlier-corrected variants)") +
  ylab("rG") +
  # ylim(0,1)+
  # gghighlight::gghighlight(label_key = Level, label_params = list(size = 8))+#, MINUTES %in% 60:90)+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  ) +
  # stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication() +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  geom_vline(xintercept = 0, size = 1, colour = "#333333") +
  # theme(title = element_text(face = 'bold'),
  #       #plot.title = element_text(color="black", size=5, face="bold"),
  #       axis.title.x = element_text(color="black", size=22),
  #       axis.title.y = element_text(color="black", size=22),
  #       legend.text = element_text(color = 'black', size = 10),
  #       legend.position = 'top',
  #       legend.title = element_text(color = 'black', size = 10))+
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  ) +
  theme(aspect.ratio = 1 / 1)


# wilcox test -------------------------------------------------------------

wilcox.test(abs(Default_df$cor_a), abs(Default_df$Between), alternative = "greater")
wilcox.test(abs(LD_df$cor_a), abs(LD_df$Between), alternative = "greater")
wilcox.test(abs(outlier_corrected_df$cor_a), abs(outlier_corrected_df$Between), alternative = "greater")

median(abs(Default_df$Between))
median(abs(Default_df$cor_a))
median(abs(Default_df$cor_a))
median(abs(outlier_corrected_df$cor_a))


# plot density ------------------------------------------------------------

ks.test(Default_df$cor_a, LD_df$cor_a, alternative = "two.sided")

p1 <- bind_cols(Default = Default, "LD corrected" = LD, "Outlier corrected" = outlier_corrected) %>%
  mutate(num = row_number()) %>%
  pivot_longer(-num, names_to = "type", values_to = "Cor_effect_size") %>%
  ggplot(aes(Cor_effect_size, color = type, fill = type)) +
  geom_density(alpha = .2, size = 2) +
  scale_fill_manual(values = c("#809F8C", "#637095", "#AB6B6B")) +
  scale_color_manual(values = c("#809F8C", "#637095", "#AB6B6B")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y = element_blank()
  ) +
  xlab("Effect Size Correlation") +
  ylab("Density") +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  geom_vline(xintercept = 0, size = 1, colour = "#333333") +
  theme_Publication() +
  theme(aspect.ratio = 1 / 1) +
  theme(legend.title = element_blank())


Correlations <-
  Correlations %>%
  mutate(horizontal = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)) %>%
    filter(outlier == TRUE)))

Correlations <-
  Correlations %>%
  mutate(b = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ))) %>%
  mutate(cor_a_all = map(b, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1])) %>%
  select(-b) %>%
  unnest(cor_a_all)


df <- Correlations %>%
  select(Phen1, Phen2, Between, cor_a, cor_a_all, horizontal)
mean <- df %>%
  mutate(diff = cor_a - cor_a_all) %>%
  summarize(mean = mean(diff), sd = sd(diff)) %>%
  pull(mean)


sd <- df %>%
  mutate(diff = cor_a - cor_a_all) %>%
  summarize(mean = mean(diff), sd = sd(diff)) %>%
  pull(sd)
df <- df %>%
  mutate(diff = cor_a - cor_a_all) %>%
  mutate(dev = ifelse(diff > mean + 3 * sd | diff < mean - 3 * sd, TRUE, FALSE))

p2 <- df %>%
  ggplot(aes(cor_a_all, cor_a, color = dev)) +
  geom_point() +
  theme_minimal() +
  # stat_function(fun = square,
  #              n = 200, color = 'steelblue4')+
  # stat_function(fun = linear,
  #               n = 200, color = 'black', size = 1.5)+
  # scale_colour_gradient2()+
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2')+
  # stat_cor(method = "pearson")
  xlab("Effect Size Correlation (Default)") +
  ylab("Effect Size Correlation (Outlier corrected)") +
  # ylim(0,1)+
  # gghighlight::gghighlight(label_key = Level, label_params = list(size = 8))+#, MINUTES %in% 60:90)+


  # stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  # bbc_style() +
  theme_Publication() +
  theme(
    axis.title.x = element_text(colour = "#809F8C", face = "bold"),
    axis.title.y = element_text(colour = "#AB6B6B", face = "bold")
  ) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  geom_vline(xintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "none") +
  annotate("label", label = "Horizontal trait pair", x = 0.57, y = -0.45, size = 4, color = "#333333") +
  annotate(
    geom = "curve", x = 0.8, y = -0.4, xend = 0.69, yend = 0.14, colour = "#333333", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    geom = "curve", x = 0.8, y = -0.4, xend = -0.5, yend = -0.05, colour = "#333333", size = 1.2,
    curvature = 0.1, arrow = arrow(length = unit(2, "mm"))
  ) +
  theme(aspect.ratio = 1 / 1) +
  scale_color_manual(
    values = c("grey40", "#e7cd79"),
    aesthetics = c("colour", "fill")
  )


# table showing horizontal trait pair and prominent horizontal loci --------


table <- df %>%
  filter(dev == "TRUE") %>%
  mutate(prominent_loci = map(horizontal, ~ .x %>%
    arrange(-abs(cor_element)) %>%
    select(num) %>%
    slice(1) %>%
    pull())) %>%
  unnest(prominent_loci) %>%
  select(Phen1, Phen2, Between, cor_a, cor_a_all, prominent_loci) # %>%
# rename(Phenotype1 = Phen1, Phenotype2 = Phen2, rG = Between, rD_default = cor_a_all, rD_outlier_corrected = cor_a, Biggest_outlier_loci = prominent_loci)


#
# CI <- table %>%
#   mutate(bootstrap = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
#                                          a2 = additive_effect %>% select(.y) %>% pull(),
#                                          num = locus_name) %>%
#                     mutate(cor_element = (a1 - mean(a1))*(a2 - mean(a2))/ sqrt(var(a1)*var(a2))) %>%
#                     mutate(H = 1.5 * IQR(cor_element)) %>%
#                     mutate(qnt1 = quantile(cor_element, probs=c(.25, .75))[1],qnt2 = quantile(cor_element, probs=c(.25, .75))[2]) %>%
#                     mutate(outlier = ifelse(cor_element >qnt2 + H | cor_element < qnt1 - H, TRUE,FALSE))
#   )) %>%
#   mutate(output =  map(bootstrap,~ bootstrap(.x, 1000) )) %>%
#   #unnest(output)
#   select(Phen1,Phen2, output) %>%
#   unnest(output) %>%
#   group_by(Phen1,Phen2) %>%
#   summarize_at(vars(cor_a), funs(!!!p_funs))
#


SE <- table %>%
  mutate(bootstrap = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  mutate(output = map(bootstrap, ~ bootstrap(.x, 1000) %>% ungroup())) %>%
  # unnest(output)
  select(Phen1, Phen2, output) %>%
  unnest(output) %>%
  group_by(Phen1, Phen2) %>%
  summarize(sd = sd(cor_a))


SE_all <- table %>%
  mutate(bootstrap = map2(Phen1, Phen2, ~ tibble(
    a1 = additive_effect %>% select(.x) %>% pull(),
    a2 = additive_effect %>% select(.y) %>% pull(),
    num = locus_name
  ) %>%
    mutate(cor_element = (a1 - mean(a1)) * (a2 - mean(a2)) / sqrt(var(a1) * var(a2))) %>%
    mutate(H = 1.5 * IQR(cor_element)) %>%
    mutate(qnt1 = quantile(cor_element, probs = c(.25, .75))[1], qnt2 = quantile(cor_element, probs = c(.25, .75))[2]) %>%
    mutate(outlier = ifelse(cor_element > qnt2 + H | cor_element < qnt1 - H, TRUE, FALSE)))) %>%
  mutate(output = map(bootstrap, ~ bootstrap_with_outlier(.x, 1000) %>% ungroup())) %>%
  # unnest(output)
  select(Phen1, Phen2, output) %>%
  unnest(output) %>%
  group_by(Phen1, Phen2) %>%
  summarize(sd = sd(cor_a))


# table %>%
#   left_join(CI, by = c('Phen1' = 'Phen1', 'Phen2' = 'Phen2')) %>%
#   select(-`50%`) %>%
#   mutate()
#   unite("CI",`2.5%`:`97.5%`, sep = ',')
#


table <- table %>%
  left_join(SE, by = c("Phen1" = "Phen1", "Phen2" = "Phen2")) %>%
  mutate(SE = round(sd * 1.96, digits = 2)) %>%
  mutate(cor_a = round(cor_a, digits = 2)) %>%
  mutate(dummy = cor_a) %>%
  unite("cor_a", dummy:SE, sep = "\u00B1") %>%
  select(-sd) %>%
  left_join(SE_all, by = c("Phen1" = "Phen1", "Phen2" = "Phen2")) %>%
  mutate(SE_all = round(sd * 1.96, digits = 2)) %>%
  mutate(cor_a_all = round(cor_a_all, digits = 2)) %>%
  mutate(dummy = cor_a_all) %>%
  unite("cor_a_all", dummy:SE_all, sep = "\u00B1") %>%
  select(-sd) %>%
  mutate(Between = round(Between, digits = 2)) %>%
  rename(
    Phenotype1 = Phen1, Phenotype2 = Phen2, rG = Between,
    rD_default = cor_a_all, rD_outlier_corrected = cor_a,
    Biggest_outlier_loci = prominent_loci
  )


write_csv(table, file = "Horizontal_trait.csv")


patch1 <- p1 + p2
patch2 <- p3 + p4 + p5

patch1 / patch2
ggsave(patch1, filename = "Fig5_1.pdf", height = 8, width = 10)
ggsave(patch2, filename = "Fig5_2.pdf", height = 8, width = 10)
