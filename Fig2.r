
library(tidyverse)
library(ggpointdensity)
library(ggpubr)
library(bbplot)
library(PhenotypeSimulator)
library(qtl)
#install.packages("MASS")                  
library("MASS")   
setwd("~/Dropbox (MIT)/Work/Research/Constraints")
source("toolbox.r")




# Fig2 --------------------------------------------------------------------

#distribution of effect size

NrSNP <- 226
rau <- runif(1,-1,1)
mu <- c(0, 0)                                   # Specify the means of the variables
sigma <- matrix(c(1, rau, rau, 1),                 # Specify the covariance matrix of the variables
                ncol = 2)

effect_1 <- draw.multivariate.laplace(NrSNP,2, 1, mu, sigma)
effect_2 <- draw.multivariate.laplace(NrSNP,2, 1, mu, sigma)

n = 20
set.seed(1010)
patch0 <- tibble(gamma = c(rep(0.5,n),rep(1,n))) %>% 
  rowwise() %>% 
  mutate(rau = runif(1,-1,1)) %>% 
  mutate(value = map2(gamma, rau, ~ draw.multivariate.laplace(NrSNP,2, .x, mu, matrix(c(1, .y, .y, 1),             
                                                                             ncol = 2))[,1])) %>% 
  ungroup() %>% 
  mutate(num = row_number()) %>% 
  unnest(a = value) %>% 
  mutate(gamma = factor(gamma)) %>% dplyr::select(gamma, a) %>% 
  group_by(gamma) %>% 
  mutate(a = scale(a)) %>% 
  #add_row(tibble(gamma = 'real', a = additive_effect$C103_A)) %>% 
  ggplot(aes(a, group = fct_rev(gamma),fill = fct_rev(gamma)))+
  geom_histogram(binwidth=.3, alpha=.5, position="identity")+
  theme_Publication()+
  theme(legend.position = 'top',
        axis.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.title =  element_text(size = 12, face = 'bold'))+
  #labs(fill = "Kurtosis")+
  xlab('Effect size')+
  ylab('Count')+
  scale_fill_manual(name = "Kurtosis", values = c("#798E87","#9C964A"), labels = c('Low', 'High'))+
  theme(aspect.ratio = 1/1)
  #scale_color_manual(values = c('#638095','#809F8C'))


high <- "#9C964A"
low <- "#798E87"

#using actual effect size 


NrSNP <- 226
sample_size <- 500

genotypes <- simulateGenotypes(N = sample_size, NrSNP = NrSNP,
                               frequencies = 0.3, # c(0.05,0.1,0.2,0.5)
                               verbose = FALSE) 
SNPinfo <- genotypes$genotypes





cross<-read.cross(format='csv',file='geno_with_pheno.csv',genotypes=c('A','B'))
class(cross)[1]<-'dh'

additive_effect <- get_additive_effect(cross)
load("Correlations.Rfile")

df <- Correlations %>% as_tibble() %>% #sample_n(10) %>%
  #mutate(Within = as.numeric(Within), Between = as.numeric(Between)) %>%
  dplyr::select(Phen1, Phen2) %>% 
  mutate(a = map2(Phen1, Phen2, ~ tibble(a1 = additive_effect %>% select(.x) %>% pull(),
                                       a2 = additive_effect %>% select(.y) %>% pull()))) %>% 
  mutate(cor_a = map(a, ~ cor(tibble(.x) %>% pull(a1), tibble(.x) %>% pull(a2))[1] )) %>%
  unnest(cor_a) %>% 
  rowwise() %>% 
  mutate(cor_g = cor(SNPinfo %*% cbind(a$a1,a$a2))[2,1] ) 


df %>%
  ggplot(aes(cor_a,cor_g))+
  geom_point(alpha = .65, color = 'Black', alpha = .1)+
  xlab('rD')+
  ylab('rG')+
  theme_Publication()+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") 


ggsave('Fig3_s1.pdf')






#polygenic
NrSNP <- 226
sample_size <- 500

genotypes <- simulateGenotypes(N = sample_size, NrSNP = NrSNP,
                               frequencies = 0.3, # c(0.05,0.1,0.2,0.5)
                               verbose = FALSE) 
SNPinfo <- genotypes$genotypes

df <- c()
for (i in c(1:1000)) {
  df <- rbind(df, simulation_for_a_given_population_exponential(
    NrSNP = NrSNP, 
    SNPinfo = SNPinfo,
    rau <- NULL,
    effect_size = FALSE, mode = 'polygenic', gamma = 1))
}

p1 <- df %>% as_tibble() %>%
  rename(genetic_cor = V1, developmental_cor = V2) %>% 
  ggplot(aes(developmental_cor,genetic_cor))+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  geom_point(alpha = .5, color = low)+
  xlab('rD')+
  ylab('rG')+
  theme_Publication()+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)
# theme(
#   #plot.title = element_text(color="black", size=5, face="bold"),
#   axis.title.x = element_text(color="black", size=30),
#   axis.title.y = element_text(color="black", size=30),
#   legend.text = element_text(color = 'black', size = 30),
#   legend.position = 'top',
#   legend.title = element_text(color = 'black', size = 12))
# annotate("label", label = "y = x", x = 0.6, y = 0.18, size = 5) +
# annotate(
#   geom = "curve", x = 0.6, y = 0.23, xend = 0.45, yend = 0.45, 
#   curvature = 0.6, arrow = arrow(length = unit(2, "mm"))
# )



#f16c23
df <- c()
for (i in c(1:1000)) {
  df <- rbind(df, simulation_for_a_given_population_exponential(
    NrSNP = NrSNP, 
    SNPinfo = SNPinfo,
    rau <- NULL,
    effect_size = FALSE, mode = 'polygenic', gamma = 0.5))
}

p2 <- df %>% as_tibble() %>%
  rename(genetic_cor = V1, developmental_cor = V2) %>% 
  ggplot(aes(developmental_cor,genetic_cor))+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  geom_point(alpha = .65, color = high)+
  xlab('rD')+
  ylab('rG')+
  theme_Publication()+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)







df <- c()
for (i in c(1:1000)) {
  df <- rbind(df, simulation_for_a_given_population_exponential(
    NrSNP = NrSNP, 
    SNPinfo = SNPinfo,
    rau <- NULL,
    effect_size = FALSE, mode = 'major', gamma = 1))
}

p3 <- df %>% as_tibble() %>%
  rename(genetic_cor = V1, developmental_cor = V2) %>% 
  ggplot(aes(developmental_cor,genetic_cor))+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  geom_point(alpha = .5, color = low)+
  xlab('rD')+
  ylab('rG')+
  theme_Publication()+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)




df <- c()
for (i in c(1:1000)) {
  df <- rbind(df, simulation_for_a_given_population_exponential(
    NrSNP = NrSNP, 
    SNPinfo = SNPinfo,
    rau <- NULL,
    effect_size = FALSE, mode = 'major', gamma = 0.5))
}

p4 <- df %>% as_tibble() %>%
  rename(genetic_cor = V1, developmental_cor = V2) %>% 
  ggplot(aes(developmental_cor,genetic_cor))+
  geom_hline(yintercept = 0, size = 1, colour="#333333") +
  geom_vline(xintercept = 0, size = 1, colour="#333333") +
  geom_point(alpha = .65, color = high)+
  xlab('rD')+
  ylab('rG')+
  theme_Publication()+
  theme(
    axis.text.y=element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y=element_blank())+
  theme(aspect.ratio = 1/1)






# 
# 
# 
# load('genotype_yeast.rdata')
# LD_threshold <- 0.5
# dummy <- t(X) %>% as_tibble() %>% 
#   mutate(num = colnames(X)) %>% 
#   rename(chr = V1) %>% 
#   nest(data = -c(chr,num)) %>% 
#   mutate(data = map(data, ~ t(.) %>% as_tibble())) 
# 
# 
# LD_cor <- expand_grid(x = dummy$num,y = dummy$num) %>%
#   left_join(dummy, by = c('x' = 'num')) %>%
#   left_join(dummy, by = c('y' = 'num')) %>%
#   filter(chr.x == chr.y & x != y) %>%
#   filter(!duplicated(paste0(pmax(x, y), pmin(x, y)))) %>%
#   mutate(LD = map2_dbl(data.x, data.y, ~ cor(.x,.y,use = "complete.obs"))) %>%
#   dplyr::select(-data.x, -data.y)
# LD <- LD_cor %>% filter(LD > LD_threshold)
# 
# 
# 
# genotype <- X[-1,]
# 
# 
# 
# df <- c()
# coef <- NULL
# 
# for (i in c(1:1000)) {
# 
#   df <- rbind(df, simulation_for_real_genotype(NrSNP = 226, rau = coef,
#                                                SNPinfo = genotype,  magnitude = 5, mode = 'polygenic', 
#                                                effect_size = FALSE, LD_block = LD, gamma = 1.0))
# }
# 
# 
# 
# p5 <- df %>% as_tibble() %>% 
#   rename(genetic_cor = V1, developmental_cor = V2) %>% 
#   ggplot(aes(developmental_cor,genetic_cor))+
#   geom_hline(yintercept = 0, size = 1, colour="#333333") +
#   geom_vline(xintercept = 0, size = 1, colour="#333333") +
#   geom_point(alpha = .5, color = '#638095')+
#   xlab('rD')+
#   ylab('rG')+
#   # xlab('Effect Size Coherence')+
#   # ylab('LD Pruned Effect Size Coherence')+
#   #theme_pubclean()+
#   theme(
#     axis.text.y=element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.ticks.y=element_blank())+
# 
#   #stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
#   #bbc_style() + 
#   theme_Publication()+
#   theme(aspect.ratio = 1/1)
# 
# 
# 
# 
# df <- c()
# coef <- NULL
# 
# 
# for (i in c(1:1000)) {
#   
#   df <- rbind(df, simulation_for_real_genotype(NrSNP = 226, rau = coef,
#                                                SNPinfo = genotype,  magnitude = 5, mode = 'polygenic', 
#                                                effect_size = FALSE, LD_block = LD, gamma = 0.5))
# }
# 
# 
# 
# p6 <- df %>% as_tibble() %>% 
#   rename(genetic_cor = V1, developmental_cor = V2) %>% 
#   ggplot(aes(developmental_cor,genetic_cor))+
#   geom_hline(yintercept = 0, size = 1, colour="#333333") +
#   geom_vline(xintercept = 0, size = 1, colour="#333333") +
#   geom_point(alpha = .65, color = '#f16c23')+
#   xlab('rD')+
#   ylab('rG')+
#   # xlab('Effect Size Coherence')+
#   # ylab('LD Pruned Effect Size Coherence')+
#   #theme_pubclean()+
#   theme(
#     axis.text.y=element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.ticks.y=element_blank())+
#   #stat_cor(method = "pearson", aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -1.0, label.y = 0.9,size = 4.5)+
#   #bbc_style() + 
#   theme_Publication()+
#   theme(aspect.ratio = 1/1)
# 

library(patchwork)
patch1 <-   p1 + p2 + grid::textGrob('No HP') + plot_layout(widths = c(2,2,1))
patch2 <-   p3 + p4 + grid::textGrob('HP') + plot_layout(widths = c(2,2,1))
#patch3 <-   p5 + p6 + grid::textGrob('LD') + plot_layout(widths = c(2,2,1))


Fig2_1 <- patch1 / patch2 
ggsave(Fig2_1, filename = 'Fig2_1.pdf')
ggsave(patch0, filename = 'Fig2_2.pdf',width = 3,height = 3)
