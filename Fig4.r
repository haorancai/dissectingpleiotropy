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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("toolbox.r")

# load data -------------------------------------------------------------------

cross<-read.cross(format='csv',file='geno_with_pheno.csv',genotypes=c('A','B'))
class(cross)[1]<-'dh'
load('LD_pruning.rdata')


# get_additive_effect -----------------------------------------------------

additive_effect <- get_additive_effect(cross)


# Default ----------------------------------------------------------------

load("Correlations.Rfile")

Correlations <- 
  Correlations  %>% 
  as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within,Between) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                         a2 = additive_effect %>% select(.y) %>% pull()))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1] )) %>%
  unnest(cor_a)

Default <- Correlations$cor_a
Default_df <- Correlations %>% 
  as_tibble() %>% 
  mutate( Between = as.numeric(Between))

p3 <- Default_df %>%
  ggplot(aes(cor_a,Between))+
  geom_point(alpha = .4, color = '#809F8C')+
  stat_function(fun = linear,
                n = 200, color = 'black', size = 1.5)+
  xlab('Effect Size Correlation \n (All variants)')+
  ylab('rG')+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  #stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication()+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  )+
  theme(aspect.ratio = 1/1)


# LD ----------------------------------------------------------------

load("Correlations.Rfile")
load('genotype_yeast.rdata')
locus_name <- colnames(X)

Correlations <-
  Correlations %>% as_tibble() %>%
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within,Between) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                         a2 = additive_effect %>% select(.y) %>% pull(),
                                         num = locus_name) %>%
                    filter(!(num %in% out_LD)))) %>%
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1] )) %>%
  unnest(cor_a)

LD<- Correlations$cor_a
  
LD_df <- Correlations %>% 
  as_tibble() %>% 
  mutate( Between = as.numeric(Between))

p4 <- LD_df %>% 
  ggplot(aes(cor_a,Between))+
  geom_point(alpha = .4, color = '#637095')+
  stat_function(fun = linear,
                n = 200, color = 'black', size = 1.5)+

  xlab('Effect Size Correlation \n (LD pruned variants)')+
  ylab('rG')+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  #stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication()+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  )+
  theme(aspect.ratio = 1/1)


# Outlier ----------------------------------------------------------------

load("Correlations.Rfile")



Correlations <- Correlations %>% as_tibble() %>% 
  mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  select(Phen1, Phen2, Within,Between) %>%
  mutate(a = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                         a2 = additive_effect %>% select(.y) %>% pull(),
                                         num = locus_name) %>%
                    mutate(cor_element = (a1 - mean(a1))*(a2 - mean(a2))/ sqrt(var(a1)*var(a2))) %>% 
                    mutate(H = 1.5 * IQR(cor_element)) %>% 
                    mutate(qnt1 = quantile(cor_element, probs=c(.25, .75))[1],qnt2 = quantile(cor_element, probs=c(.25, .75))[2]) %>% 
                    mutate(outlier = ifelse(cor_element >qnt2 + H | cor_element < qnt1 - H, TRUE,FALSE)) %>% 
                    filter(outlier == FALSE)
  )) %>% 
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1] )) %>%
  unnest(cor_a)

outlier_corrected <- Correlations$cor_a
outlier_corrected_df <- Correlations %>% 
  as_tibble() %>% 
  mutate(Between = as.numeric(Between))

p5 <- outlier_corrected_df %>% 
  ggplot(aes(cor_a,Between))+
  geom_point(alpha = .4, color = '#AB6B6B')+
  stat_function(fun = linear,
                n = 200, color = 'black', size = 1.5)+
  xlab('Effect Size Correlation \n (Outlier-corrected variants)')+
  ylab('rG')+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  #stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
  theme_Publication()+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  annotate("label", label = "y = x", x = 0.8, y = 0.18, size = 5) +
  annotate(
    geom = "curve", x = 0.8, y = 0.23, xend = 0.50, yend = 0.45, colour = "black", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  )+
  theme(aspect.ratio = 1/1)






# wilcox test -------------------------------------------------------------

wilcox.test(abs(Default_df$cor_a),abs(Default_df$Between), alternative = 'greater')
wilcox.test(abs(LD_df$cor_a),abs(LD_df$Between), alternative = 'greater')
wilcox.test(abs(outlier_corrected_df$cor_a),abs(outlier_corrected_df$Between), alternative = 'greater')

median(abs(Default_df$Between))
median(abs(Default_df$cor_a))
median(abs(Default_df$cor_a))
median(abs(outlier_corrected_df$cor_a))


# plot density ------------------------------------------------------------

ks.test(Default_df$cor_a, LD_df$cor_a, alternative = 'two.sided')

p1 <- bind_cols(Default = Default, 'LD corrected'= LD, 'Outlier corrected'= outlier_corrected) %>% 
  mutate(num = row_number()) %>% 
  pivot_longer(-num, names_to = 'type', values_to = 'Cor_effect_size') %>% 
  ggplot(aes(Cor_effect_size, color = type, fill = type))+ geom_density(alpha = .2, size = 2)+
  scale_fill_manual(values = c('#809F8C','#637095','#AB6B6B'))+
  scale_color_manual(values = c('#809F8C','#637095','#AB6B6B'))+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  xlab('Effect Size Correlation')+
  ylab('Density')+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  theme_Publication()+
  theme(aspect.ratio = 1/1)+
  theme(legend.title = element_blank())



Correlations <- 
  Correlations %>% 
  mutate(horizontal = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                                  a2 = additive_effect %>% select(.y) %>% pull(),
                                                  num = locus_name) %>%
                             mutate(cor_element = (a1 - mean(a1))*(a2 - mean(a2))/ sqrt(var(a1)*var(a2))) %>% 
                             mutate(H = 1.5 * IQR(cor_element)) %>% 
                             mutate(qnt1 = quantile(cor_element, probs=c(.25, .75))[1],qnt2 = quantile(cor_element, probs=c(.25, .75))[2]) %>% 
                             mutate(outlier = ifelse(cor_element >qnt2 + H | cor_element < qnt1 - H, TRUE,FALSE)) %>% 
                             filter(outlier == TRUE)
  )) 

Correlations <- 
  Correlations %>% 
  mutate(b = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                         a2 = additive_effect %>% select(.y) %>% pull(),
                                         num = locus_name) 
  ))  %>% 
  mutate(cor_a_all = map(b, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1] )) %>%
  select(-b) %>% 
  unnest(cor_a_all)


df <- Correlations %>% 
  select(Phen1, Phen2,Between, cor_a,cor_a_all, horizontal)
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
  mutate(dev = ifelse(diff > mean + 3*sd | diff < mean - 3*sd, TRUE, FALSE)) 

p2 <- df %>% 
  ggplot(aes(cor_a_all,cor_a, color = dev))+geom_point()+theme_minimal()+
  xlab('Effect Size Correlation (Default)')+
  ylab('Effect Size Correlation (Outlier corrected)')+
  theme_Publication()+
  theme(axis.title.x = element_text(colour = "#809F8C", face = 'bold'),
        axis.title.y = element_text(colour = "#AB6B6B", face = 'bold'))+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  theme(legend.position = 'none')+
  annotate("label", label = "Horizontal trait pair", x = 0.57, y = -0.45, size = 4, color = "#333333") +
  annotate(
    geom = "curve", x = 0.8, y = -0.4, xend = 0.69, yend = 0.14, colour = "#333333", size = 1.2,
    curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
  )+
  annotate(
    geom = "curve", x = 0.8, y = -0.4, xend = -0.5, yend = -0.05, colour = "#333333", size = 1.2,
    curvature = 0.1, arrow = arrow(length = unit(2, "mm"))
  )+
  theme(aspect.ratio = 1/1)+
  scale_color_manual(values = c("grey40","#e7cd79"),
                     aesthetics = c("colour", "fill"))





# table showing horizontal trait pair and prominent horizontal loci --------



table<- df %>% filter(dev == 'TRUE') %>% 
  mutate(prominent_loci = map(horizontal, ~ .x %>% arrange(-abs(cor_element)) %>% select(num) %>% slice(1) %>% pull())) %>% 
  unnest(prominent_loci) %>% 
  select(Phen1,Phen2, Between, cor_a,cor_a_all, prominent_loci) %>% 
  rename(Phenotype1 = Phen1, Phenotype2 = Phen2, rG = Between, rD_default = cor_a_all, rD_outlier_corrected = cor_a, Biggest_outlier_loci = prominent_loci)


write_csv(table, file = 'Horizontal_trait.csv')


patch1 <-  p1 +  p2
patch2 <- p3 + p4 + p5

patch1 / patch2
ggsave(patch1, filename = 'Fig4_1.pdf', height = 8, width = 10)
ggsave(patch2, filename = 'Fig4_2.pdf', height = 8, width = 10)



