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
library(car)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("function.r")
# load data -------------------------------------------------------------------


load("merged_correlation.Rfile")
load('genotype_yeast.rdata')
load("Correlations.Rfile")
cross<-read.cross(format='csv',file='geno_with_pheno.csv',genotypes=c('A','B'))
class(cross)[1]<-'dh'


# get_additive_effect -----------------------------------------------------

additive_effect <- get_additive_effect(cross)
locus_name <- colnames(X)

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


Merged_Correlations<- Merged_Correlations %>% 
  mutate(Percent = case_when(Percent == '0uM_GdA' ~ '0µM_GdA',
                             Percent == '8.5uM_GdA' ~ '8.5µM_GdA',
                             Percent == '25uM_GdA' ~ '25µM_GdA',
                             Percent =='100uM_GdA' ~ '100µM_GdA'))
  
  
Merged_Correlations <- Merged_Correlations %>% select(-Between, -Within) %>% 
  left_join(Correlations %>% select(Phen1, Phen2, cor_a,Between), 
            by = c('Phen1' = 'Phen1', 'Phen2' = 'Phen2')) %>% 
  rename(rG = Between, change = B_Difference, rD = cor_a)



#accounting for 0uM GdA 
Merged_Correlations %>% as_tibble() %>% #sample_n(10) %>%
  mutate(rG = as.numeric(rG)) %>%
  filter(Percent != '0µM_GdA') %>% 
  left_join(Merged_Correlations %>% filter(Percent == '0µM_GdA') %>% select(Pair, change), by = c('Pair' = 'Pair')) %>% 
  mutate(change = abs(change.x - change.y)) %>% 
  ggplot(aes(abs(rD), change, color = abs(rG)))+
  ylab('changes of rG')+
  xlab('rD')+
  geom_point()+
  scale_color_viridis(name = 'rG',option="magma") +
  facet_grid(cols = vars(Percent))+
  theme_Publication()+
  theme(legend.text = element_text(size = 10), legend.key.width = unit(2, 'cm'))+
  theme(aspect.ratio = 1/1)

ggsave(filename = 'Fig5_2.png')


data <- Merged_Correlations %>% as_tibble() %>% #sample_n(10) %>%
  mutate(rG = as.numeric(rG)) %>%
  filter(Percent != '0µM_GdA') %>% 
  left_join(Merged_Correlations %>% filter(Percent == '0µM_GdA') %>% select(Pair, change), by = c('Pair' = 'Pair')) %>% 
  mutate(change = change.x - change.y) %>% 
  select(change,rD,rG, Percent)



# regression results

fit <- data %>% 
  as_tibble() %>% 
  #filter(Percent == "25uM_GdA") %>% 
  select(change,rD,rG) %>% 
  mutate_all(abs) %>% 
  {lm(change ~ rD + rG, data = .)} 


#plot_summs(fit1, robust = TRUE, inner_ci_level = .9)

fit1 <- data %>% 
  as_tibble() %>% 
  filter(Percent == "8.5µM_GdA") %>% 
  select(change,rD,rG) %>% 
  mutate_all(abs) %>% 
  {lm(change ~ rD + rG, data = .)} 



fit2 <- data %>% 
  as_tibble() %>% 
  filter(Percent == "25µM_GdA") %>% 
  select(change,rD,rG) %>% 
  mutate_all(abs) %>% 
  {lm(change ~ rD + rG, data = .)} 



fit3 <- data %>% 
  as_tibble() %>% 
  filter(Percent == "100µM_GdA") %>% 
  select(change,rD,rG) %>% 
  mutate_all(abs) %>% 
  {lm(change ~ rD + rG, data = .)} 

vif(fit2)
library(broom.mixed)
p0 <- plot_summs(fit,fit1,fit2,fit3,
                 model.names = c("Aggregated", "8.5µM GdA","25µM GdA","100µM GdA"))
#ggsave(filename = 'mlm_summ.pdf')



  


# bivariate normal sampling

rbvn<-function (n, X, rho)
{
  
  X1 <- rnorm(n, mean(X) +  rho *
                (X - mean(X)), sqrt((1 - rho^2)*sd(X)^2))
  X1
}

set.seed(1010)
null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "8.5µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        #rowwise() %>% 
                        mutate(rGn = rbvn(nrow(data)/3,rG,0.93)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p1 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "8.5µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05, color = 'black')+
  geom_point(alpha = .7) +
  ylab('Slope Estimate (rG)')+
  xlab('Slope Estimate (rD)')+
  labs(title="8.5µM GdA")+
  # xlim(-0.5, 0.5)+
  # ylim(-0.5, 0.25)+
  scale_color_manual(values = c('#D9091E','black'))+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p1 <- ggExtra::ggMarginal(p1, margin = 'x', groupColour = TRUE, groupFill = TRUE)


set.seed(1010)
null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "25µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        #rowwise() %>% 
                        mutate(rGn = rbvn(nrow(data)/3,rG,0.93)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p2 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "25µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05, color = 'black')+
  geom_point(alpha = .7) +
  ylab('Slope Estimate (rG)')+
  xlab('Slope Estimate (rD)')+
  labs(title="25µM GdA")+
  # xlim(-0.5, 0.5)+
  # ylim(-0.5, 0.25)+
  scale_color_manual(values = c('#D9091E','black'))+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p2 <- ggExtra::ggMarginal(p2, margin = 'x', groupColour = TRUE, groupFill = TRUE)


set.seed(1010)
null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "100µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        #rowwise() %>% 
                        mutate(rGn = rbvn(nrow(data)/3,rG,0.93)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p3 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "100µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05, color = 'black')+
  scale_color_manual(values = c('#D9091E','black'))+
  geom_point(alpha = .6) +
  ylab('Slope Estimate (rG)')+
  xlab('Slope Estimate (rD)')+
  labs(title="100µM GdA")+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p3 <- ggExtra::ggMarginal(p3, margin = 'x', groupColour = TRUE, groupFill = TRUE)




rG_rGn_bvn <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)




cowplot::plot_grid(p0+theme(legend.position=c(0.75,0.9), 
                            legend.title = element_blank(),
                            legend.text = element_text(size = 10, face = 'bold')), 
                   rG_rGn_bvn,align = "h",  axis = "b", rel_widths = c(1, 2) , ncol = 2, nrow = 1)#,labels=c("b", "c")
cowplot::ggsave2(filename = 'Fig5_1.pdf', height = 5, width = 16)


# SI ----------------------------------------------------------------------


null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "8.5µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        rowwise() %>% 
                        mutate(rGn = rG + rnorm(1,mean = 0, sd = 0.1)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p1 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "8.5µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05,color = 'black')+
  geom_point(alpha = .7) +
  ylab('Estimate (rG)')+
  xlab('Estimate (rD)')+
  labs(title="8.5µM GdA")+
  # xlim(-0.5, 0.5)+
  # ylim(-0.5, 0.25)+
  scale_color_manual(values = c('#D9091E','black'))+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p1 <- ggExtra::ggMarginal(p1, margin = 'x', groupColour = TRUE, groupFill = TRUE)

null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "25µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        rowwise() %>% 
                        mutate(rGn = rG + rnorm(1,mean = 0, sd = 0.1)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p2 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "25µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05, color = 'black')+
  geom_point(alpha = .7) +
  ylab('Estimate (rG)')+
  xlab('Estimate (rD)')+
  labs(title="25µM GdA")+
  scale_color_manual(values = c('#D9091E','black'))+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p2 <- ggExtra::ggMarginal(p2, margin = 'x', groupColour = TRUE, groupFill = TRUE)


null <- tibble(num = c(1:200)) %>% 
  mutate(output = map(num, ~ data %>% 
                        as_tibble() %>% 
                        filter(Percent == "100µM_GdA") %>% 
                        select(change,rD,rG) %>% 
                        rowwise() %>% 
                        mutate(rGn = rG + rnorm(1,mean = 0, sd = 0.1)) %>% 
                        mutate_all(abs) %>% 
                        {lm(change ~ rGn + rG, data = .)} %>% 
                        tidy() %>% 
                        select(term, estimate, std.error) %>% 
                        pivot_longer(cols = -1) %>% 
                        pivot_wider(names_from = 'term',values_from = "value") %>% 
                        select(name, rGn, rG) %>% 
                        pivot_wider(names_from = 'name', values_from = c('rGn', 'rG')))
  )


p3 <- null %>% unnest(output) %>% select(-num) %>% 
  add_row(data %>% 
            as_tibble() %>% 
            filter(Percent == "100µM_GdA") %>% 
            select(change,rD,rG)  %>% 
            mutate_all(abs) %>% 
            {lm(change ~ rD + rG, data = .)} %>% 
            tidy() %>% 
            select(term, estimate, std.error) %>% 
            pivot_longer(cols = -1) %>% 
            pivot_wider(names_from = 'term',values_from = "value") %>% 
            select(name, rD, rG) %>% 
            rename(rGn = rD) %>% 
            pivot_wider(names_from = 'name', values_from = c('rGn', 'rG'))) %>% 
  mutate(xmin = rGn_estimate - 2* rGn_std.error, xmax = rGn_estimate + 2* rGn_std.error,
         ymin = rG_estimate - 2* rG_std.error, ymax = rG_estimate + 2* rG_std.error) %>% 
  mutate(factor = c(rep('1',200), '0')) %>% 
  ggplot(aes(rGn_estimate,rG_estimate,color = factor))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = ymin,ymax = ymax), alpha = .05, color = 'black') + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax ), alpha = .05, color = 'black')+
  geom_point(alpha = .7) +
  labs(title="100µM GdA")+
  scale_color_manual(values = c('#D9091E','black'))+
  theme_Publication()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank())


p3 <- ggExtra::ggMarginal(p3, margin = 'x', groupColour = TRUE, groupFill = TRUE)




rG_rGn_n <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)
ggsave(filename = 'Fig5_SI.pdf', rG_rGn_n)

