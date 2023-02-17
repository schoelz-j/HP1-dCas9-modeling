# qrt-pcr data from dCas9 studies

# Analyze all datasets, one at a time, dealing withoutliers and producing
# visualizations/statistics

library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
getwd()

# Cpr11B
cpr11b_data <- read.csv('cpr11b_data.csv')
str(cpr11b_data)

# Create a sample ID column
cpr11b_data <- cpr11b_data %>%
  mutate('Sample.ID' = paste(Genotype, Condition, Bio.rep, sep = '_'))
# Check - if the data is been entered correctly, then each sample ID should
# appear 4 times.
check1 <- cpr11b_data %>%
  group_by(Sample.ID) %>%
  summarise(total = n())
# This looks good.

# I only want samples from the second qPCR run (leave out the first run)
#cpr11b_data <- cpr11b_data %>%
#filter(Sample.ID != 'GFP_Cpr11B_1') %>%
#filter(Sample.ID != 'GFP_Cpr11B_2') %>%
#filter(Sample.ID != 'GFP_Cpr11B_3') %>%
#filter(Sample.ID != 'HP1B_Cpr11B_1') %>%
#filter(Sample.ID != 'HP1B_Cpr11B_2') %>%
#filter(Sample.ID != 'HP1B_Cpr11B_3') %>%
#filter(Sample.ID != 'HP1C_Cpr11B_1') %>%
#filter(Sample.ID != 'HP1C_Cpr11B_2') %>%
#filter(Sample.ID != 'HP1C_Cpr11B_3') 


# I'm happy with how this looks, now I want to write functions to make it
# quicker to run the code on remaining datasets.
# I want one function to generate the final, analyzed table and one function
# that outputs my plot.

# table function
dd_results <- function(rt_dataframe, target_gene){
  
  # Calculate means of technical replicates
  rep_means <- rt_dataframe %>%
    group_by(Sample.ID, target.gene) %>%
    summarise(mean.Cq = mean(Cq, na.rm = T)) %>%
    data.frame()  
  
  # Calculate relative expression within each sample
  rpl32 <- rep_means %>%
    filter(target.gene == 'RpL32')
  delta1s <- left_join(rep_means, rpl32, by = 'Sample.ID') %>%
    filter(target.gene.x == target_gene) %>%
    mutate(delta1 = mean.Cq.x - mean.Cq.y) %>%
    separate(Sample.ID, c('Genotype', 'Condition', 'Bio.Rep'))
  
  # Calculate double delta values
  baselines <- delta1s %>%
    filter(Condition == "Baseline") %>%
    group_by(target.gene.x, Genotype) %>%
    summarise(base_d1 = mean(delta1, na.rm = T)) %>%
    data.frame() %>%
    mutate(gene_geno = paste(Genotype, target.gene.x, sep = '_'))
  
  tethers <- delta1s %>%
    filter(Condition != "Baseline") %>%
    mutate(gene_geno = paste(Genotype, target.gene.x, sep = '_'))
  
  delta_deltas <- left_join(tethers, baselines, by = 'gene_geno') %>%
    mutate(double_delt = delta1 - base_d1) %>%
    mutate(fold_change = 2^(-1*double_delt)) %>%
    filter(target.gene.x.x != 'RpL32') %>%
    select(Genotype.x, Bio.Rep, delta1, base_d1,
           double_delt, fold_change)
}

# Did this work?
test <- dd_results(cpr11b_data, 'Cpr11B')
#identical(test, delta_deltas)

# Plot function
rt_plot <- function(dd_dataframe, target.gene){
  
  # Prep for geom_segment
  tether_means <- dd_dataframe %>%
    group_by(Genotype.x) %>%
    summarise(average_fc = mean(fold_change))
  
  ggplot(dd_dataframe, aes(x = Genotype.x, y = fold_change, color = Genotype.x))+
    geom_jitter(size = 2, alpha = 0.25, width = 0.4)+
    stat_summary(fun = mean, geom = "point", size = 5)+
    theme_classic()+
    coord_flip()+
    geom_hline(yintercept = 1, lwd = 1)+
    xlab("Tether")+
    ylab("Fold Change")+
    scale_color_manual(values = c('grey', '#473E92', '#2D8775', '#8F3086'))+
    ggtitle(target.gene)+
    theme(legend.position = 'none')+
    geom_segment(data = tether_means, aes(x = Genotype.x, y = average_fc,
                                          xend = Genotype.x, yend = 1), lwd = 1.5)+
    scale_y_continuous(expand = c(0,0))
}

# Did this work?
test <- dd_results(cpr11b_data, 'Cpr11B')
#identical(test, delta_deltas)
cpr11b_results <- rt_plot(test, 'Cpr11B')
cpr11b_results
# Stats
summary(aov(test$double_delt ~ test$Genotype.x))
pairwise.t.test(test$double_delt, test$Genotype.x,
                p.adjust.method = 'fdr')

# Work through with remaining genes
dgt3_data <- read.csv('dgt3_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(Sample.ID != 'HP1C_Dgt3_6')
dgt3_dd <- dd_results(dgt3_data, 'Dgt3')
dgt3_results <- rt_plot(dgt3_dd, 'Dgt3')
summary(aov(dgt3_dd$double_delt ~ dgt3_dd$Genotype.x))
pairwise.t.test(dgt3_dd$double_delt, dgt3_dd$Genotype.x,
                p.adjust.method = 'fdr')

dgt3_results

cecA1_data <- read.csv('cecA1_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) #%>%
  filter(Sample.ID != 'GFP_Baseline_1') %>%
  filter(Sample.ID != 'HP1B_Baseline_3') %>%
  filter(Sample.ID != 'HP1C_Baseline_2') %>%
  filter(Sample.ID != 'GFP_CecA1_1') %>%
  filter(Sample.ID != 'GFP_CecA1_2') %>%
  filter(Sample.ID != 'HP1C_CecA1_3') %>%
  filter(Sample.ID != 'HP1C_CecA1_8') %>%
  filter(Sample.ID != 'HP1B_CecA1_6') %>%
  filter(Sample.ID != 'HP1B_CecA1_7') %>%
  filter(Sample.ID != 'HP1C_CecA1_4') %>%
  filter(Sample.ID != 'HP1C_CecA1_5') %>%
  filter(Sample.ID != "HP1C_CecA1_8") %>%
  filter(Sample.ID != 'GFP_CecA1_5') %>%
  filter(Sample.ID != 'GFP_CecA1_6') %>%
  filter(Sample.ID != 'HP1B_Baseline_10') %>%
  filter(Sample.ID != 'Baseline_GFP_11') %>%
  filter(Sample.ID != 'HP1a_CecA1_1') #%>%
  filter(Sample.ID != 'HP1a_CecA1_4') #%>%
  filter(Sample.ID != 'GFP_CecA1_3')

cecA1_clean <- read.csv('cecA1_clean.csv') %>%
  mutate(Sample.ID = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(condition == 'Baseline' |
           date == 3) %>%
  filter(Sample.ID != 'GFP_CecA1_13')
cecA1_dd <- dd_results(cecA1_clean, 'CecA1')
ceca1_results <- rt_plot(cecA1_dd, 'CecA1')
summary(aov(cecA1_dd$double_delt ~ cecA1_dd$Genotype.x))
pairwise.t.test(cecA1_dd$double_delt, cecA1_dd$Genotype.x,
                p.adjust.method = 'fdr')
ceca1_results

mtk_data <- read.csv('mtk_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(Sample.ID != 'HP1a_mtk_3') %>%
  filter(Sample.ID != 'HP1a_mtk_2')
mtk_dd <- dd_results(mtk_data, 'mtk')  
mtk_results <- rt_plot(mtk_dd, 'Mtk')
summary(aov(mtk_dd$double_delt ~ mtk_dd$Genotype.x))
pairwise.t.test(mtk_dd$double_delt, mtk_dd$Genotype.x,
                p.adjust.method = 'fdr')
mtk_results

##########################################
# pyr date 2 = rerun of samples by A.K.
# date 1 = same samples, different run by J.S.
pyr_data <- read.csv('pyr_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(date == 2)
pyr_dd <- dd_results(pyr_data, 'pyr')
rt_plot(pyr_dd, 'pyr')
t.test(pyr_dd$double_delt[pyr_dd$Genotype.x=='HP1a'],
       pyr_dd$double_delt[pyr_dd$Genotype.x=='GFP'])

# egr date 2 - qpcr run by A.K.
# removed one outlier from GFP condition
egr_data <- read.csv('egr_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(Sample.ID != 'GFP_egr_5') %>%
  filter(date == 2)
egr_dd <- dd_results(egr_data, 'egr')
rt_plot(egr_dd, 'egr')
t.test(egr_dd$double_delt[egr_dd$Genotype.x=='HP1a'],
       egr_dd$double_delt[egr_dd$Genotype.x=='GFP'])

# mats date 2 - qpcr run by A.K.
mats_data <- read.csv('mats_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  #filter(date == 2)
mats_dd <- dd_results(mats_data, 'mats')
rt_plot(mats_dd, 'mats')
t.test(mats_dd$double_delt[mats_dd$Genotype.x=='HP1a'],
       mats_dd$double_delt[mats_dd$Genotype.x=='GFP'])

################################################ 

############ HP1B tethering
shaw_data <- read.csv('shaw_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  
shaw_dd <- dd_results(shaw_data, 'Shaw')  
rt_plot(shaw_dd, 'Shaw')+
  scale_color_manual(values = c('#BEBEBE', '#2D8775'))
t.test(shaw_dd$double_delt[shaw_dd$Genotype.x=='HP1B'],
       shaw_dd$double_delt[shaw_dd$Genotype.x=='GFP'])

# A.K. reran Cg76 baseline expression for HP1B genotype, use
# those measurements for this analysis (remove old baseline measurements)
cg76_data <- read.csv('cg76-data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(Sample.ID != 'HP1B_Baseline_1') %>%
  filter(Sample.ID != 'HP1B_Baseline_2') %>%
  filter(Sample.ID != 'HP1B_Baseline_3') %>%
  filter(Sample.ID != 'HP1B_Baseline_4')
cg76_dd <- dd_results(cg76_data, 'CG13676')
rt_plot(cg76_dd, 'CG13676')+
  scale_color_manual(values = c('#BEBEBE', '#2D8775'))
t.test(cg76_dd$double_delt[cg76_dd$Genotype.x == 'HP1B'],
       cg76_dd$double_delt[cg76_dd$Genotype.x == 'GFP'])

rab3_data <- read.csv('rab3-data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
rab3_dd <- dd_results(rab3_data, 'Rab3-GEF')
rt_plot(rab3_dd, 'Rab3-GEF')+
  scale_color_manual(values = c('#BEBEBE', '#2D8775'))

t.test(rab3_dd$double_delt[rab3_dd$Genotype.x == 'HP1B'],
       rab3_dd$double_delt[rab3_dd$Genotype.x == 'GFP'])

########## HP1C tethering
# Day 2 = processed by A.K.
alpha_data <- read.csv('alphatubulin_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
  filter(date == 2)
alpha_dd <- dd_results(alpha_data, 'AlphaTubulin84b')
rt_plot(alpha_dd, 'Alpha Tubulin 84b')+
  scale_color_manual(values = c('#BEBEBE', '#8F3086'))
t.test(alpha_dd$double_delt[alpha_dd$Genotype.x == 'HP1C'],
       alpha_dd$double_delt[alpha_dd$Genotype.x == 'GFP'])
# NCR changed x == 'HP1B' to 'HP1C' to make this work 2/17/22. 

# Use new GFP baseline processed by A.K.
cg26_data <- read.csv('CG13326_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_'))
cg26_dd <- dd_results(cg26_data, 'CG13326')
rt_plot(cg26_dd, 'CG13326')+
  scale_color_manual(values = c('#BEBEBE', '#8F3086'))
t.test(cg26_dd$double_delt[cg26_dd$Genotype.x=='HP1C'],
       cg26_dd$double_delt[cg26_dd$Genotype.x=='GFP'])


# dpp - data entirely processed by A.K.
dpp_data <- read.csv('dpp_data.csv') %>%
  mutate('Sample.ID' = paste(Genotype, condition, bio.rep, sep = '_')) %>%
dpp_dd <- dd_results(dpp_data, 'dpp')
rt_plot(dpp_dd, 'Dpp')+
  scale_color_manual(values = c('#BEBEBE', '#8F3086'))

t.test(dpp_dd$double_delt[dpp_dd$Genotype.x == 'HP1C'],
       dpp_dd$double_delt[dpp_dd$Genotype.x == 'GFP'])

# Code below is an alternative analysis, normalize expression to GFP instead of
# to baseline, then plot baseline expression and post-tether expression

# New function for analysis step

gfp_dd_results <- function(rt_dataframe, target_gene){
  rep_means <- rt_dataframe %>%
    group_by(Sample.ID, target.gene) %>%
    summarise(mean.Cq = mean(Cq, na.rm = T)) %>%
    data.frame()  
  
  # Calculate relative expression within each sample
  rpl32 <- rep_means %>%
    filter(target.gene == 'RpL32')
  delta1s <- left_join(rep_means, rpl32, by = 'Sample.ID') %>%
    filter(target.gene.x == target_gene) %>%
    mutate(delta1 = mean.Cq.x - mean.Cq.y) %>%
    separate(Sample.ID, c('Genotype', 'Condition', 'Bio.Rep'))
  
  gfp <- delta1s %>%
    filter(Genotype == 'GFP') %>%
    group_by(Condition) %>%
    summarise(mean_delta = mean(delta1, na.rm =T))
  
  hp1 <- delta1s %>%
    filter(Genotype != 'GFP') %>%
    select(Genotype, Condition, Bio.Rep, delta1)
  
  
  calc_double_delt <- function(i){
    double_delt <- hp1$delta1[i] - gfp$mean_delta[gfp$Condition == hp1$Condition[i]]
    double_delt
  }
  hp1$double_delt <- unlist(lapply(1:nrow(hp1), calc_double_delt))
  hp1 <- hp1 %>%
    mutate(fold_change = 2^(-1 * double_delt))
  hp1
}

summarize_qpcr <- function(qpcr){
  qpcr %>%
    group_by(Genotype, Condition) %>%
    summarise(mean_fc = mean(fold_change, na.rm = T),
              sd_fc = sd(fold_change, na.rm = T),
              mean_dd = mean(double_delt, na.rm = T),
              sd_dd = sd(double_delt, na.rm = T))
}

cpr11b_gfp <- gfp_dd_results(cpr11b_data, 'Cpr11B')
summarize_qpcr(cpr11b_gfp)
dgt3_gfp <- gfp_dd_results(dgt3_data, 'Dgt3')
summarize_qpcr(dgt3_gfp)
cecA1_gfp <- gfp_dd_results(cecA1_data, 'CecA1')
summarize_qpcr(cecA1_gfp)
mtk_gfp <- gfp_dd_results(mtk_data, 'mtk')
summarize_qpcr(mtk_gfp)

mats_gfp <- gfp_dd_results(mats_data, 'mats')
summarize_qpcr(mats_gfp)
pyr_gfp <- gfp_dd_results(pyr_data, 'pyr')
summarize_qpcr(pyr_gfp)
egr_gfp <- gfp_dd_results(egr_data, 'egr')
summarize_qpcr(egr_gfp)

alpha_gfp <- gfp_dd_results(alpha_data, 'AlphaTubulin84b')
summarize_qpcr(alpha_gfp)
cg76_gfp <- gfp_dd_results(cg76_data, 'CG13676')
summarize_qpcr(cg76_gfp)
rab3_gfp <- gfp_dd_results(rab3_data, 'Rab3-GEF')
summarize_qpcr(rab3_gfp)

cg26_gfp <- gfp_dd_results(cg26_data, 'CG13326')
summarize_qpcr(cg26_gfp)
shaw_gfp <- gfp_dd_results(shaw_data, 'shaw')
summarize_qpcr(shaw_gfp)
dpp_gfp <- gfp_dd_results(dpp_data, 'dpp')
summarize_qpcr(dpp_gfp)
#dpp - this gene needs to be redone

# Looking at these, I think what I want to do is to do a combination,
# normalize the change in HP1 treatments to the change in GFP control

str(dgt3_dd)

gfp_mean <- dgt3_dd %>%
  filter(Genotype.x == 'GFP') %>%
  summarise(gfp_val = mean(double_delt))
gfp_mean$gfp_val
hp1 <- dgt3_dd %>%
  filter(Genotype.x != 'GFP') %>%
  mutate(GFP_Normalized = double_delt - gfp_mean$gfp_val[1]) %>%
  mutate(fold_change_gfp = 2^(-1*GFP_Normalized))
rt_plot(hp1, 'Dgt3')

gfp_norm <- function(rtqpcr_dd){
  gfp_mean <- rtqpcr_dd %>%
    filter(Genotype.x == 'GFP') %>%
    summarise(gfp_val = mean(double_delt, na.rm = T))
  
  hp1 <- rtqpcr_dd %>%
    filter(Genotype.x != 'GFP') %>%
    mutate(GFP_Normalized = double_delt - gfp_mean$gfp_val[1]) %>%
    mutate(fold_change_gfp = 2^(-1*GFP_Normalized))
}
dgt3_gfp <- gfp_norm(dgt3_dd)

gfp_plot <- function(dd_dataframe, target_gene){
  
  # Prep for geom_segment
  tether_means <- dd_dataframe %>%
    group_by(Genotype.x) %>%
    summarise(average_fc = mean(fold_change_gfp, na.rm = T))
  
  ggplot(dd_dataframe, aes(x = Genotype.x, y = fold_change_gfp, color = Genotype.x))+
    geom_jitter(size = 2, alpha = 0.25, width = 0.4)+
    stat_summary(fun = mean, geom = "point", size = 5)+
    theme_classic()+
    coord_flip()+
    geom_hline(yintercept = 1, color = 'green')+
    xlab("Tether")+
    ylab("Fold Change")+
    scale_fill_viridis_d(option = 'E')+
    ggtitle(target_gene)+
    theme(legend.position = 'none')+
    geom_segment(data = tether_means, aes(x = Genotype.x, y = average_fc,
                                          xend = Genotype.x, yend = 1))+
    scale_y_continuous(expand = c(0,0))
}
gfp_plot(dgt3_gfp, 'Dgt3')

cpr11b_dd <- dd_results(cpr11b_data, 'Cpr11B')
cpr11b_gfp <- gfp_norm(cpr11b_dd)
gfp_plot(cpr11b_gfp, 'Cpr11B')

cecA1_gfp <- gfp_norm(cecA1_dd)
gfp_plot(cecA1_gfp, 'CecA1')

mtk_gfp <- gfp_norm(mtk_dd)
gfp_plot(mtk_gfp, 'Mtk')

mats_gfp <- gfp_norm(mats_dd)
pyr_gfp <- gfp_norm(pyr_dd)
egr_gfp <- gfp_norm(egr_dd)

mats_gfp$gene = rep('mats', nrow(mats_gfp))
pyr_gfp$gene = rep('pyr', nrow(pyr_gfp))
egr_gfp$gene = rep('egr', nrow(egr_gfp))

hp1a_tethers <- data.frame(rbind(mats_gfp, pyr_gfp, egr_gfp))
str(hp1a_tethers)

group_plot <- function(dd_dataframe){
  
  # Prep for geom_segment
  tether_means <- dd_dataframe %>%
    group_by(gene) %>%
    summarise(average_fc = mean(fold_change))
  
  ggplot(dd_dataframe, aes(x = gene, y = fold_change, color = gene))+
    geom_jitter(size = 2, alpha = 0.25, width = 0.4)+
    stat_summary(fun = mean, geom = "point", size = 5)+
    theme_classic()+
    coord_flip()+
    geom_hline(yintercept = 1, color = 'green')+
    xlab("Tether")+
    ylab("Fold Change")+
    scale_fill_viridis_d(option = 'E')+
    #ggtitle(target_gene)+
    theme(legend.position = 'none')+
    geom_segment(data = tether_means, aes(x = gene, y = average_fc,
                                          xend = gene, yend = 1))+
    scale_y_continuous(expand = c(0,0))
}
group_plot(hp1a_tethers)

alpha_gfp <- gfp_norm(alpha_dd)
alpha_gfp$gene <- rep('AlphaTubulin84b', nrow(alpha_gfp))
cg76_gfp <- gfp_norm(cg76_dd)
cg76_gfp$gene <- rep('CG13676', nrow(cg76_gfp))
rab3_gfp <- gfp_norm(rab3_dd)
rab3_gfp$gene <- rep('Rab3-GEF', nrow(rab3_gfp))
str(alpha_gfp)
str(cg76_gfp)
str(rab3_gfp)
hp1b_tethers <- data.frame(rbind(
  alpha_gfp,
  cg76_gfp,
  rab3_gfp
))
group_plot(hp1b_tethers)

cg26_gfp <- gfp_norm(cg26_dd)
shaw_gfp <- gfp_norm(shaw_dd)
dpp_gfp <- gfp_norm(dpp_dd)
cg26_gfp$gene <- rep('CG13326', nrow(cg26_gfp))
shaw_gfp$gene <- rep('Shaw', nrow(shaw_gfp))
dpp_gfp$gene <- rep('Dpp', nrow(dpp_gfp))

hp1c_tethers <- data.frame(rbind(
  cg26_gfp,
  shaw_gfp,
  dpp_gfp
))
group_plot(hp1c_tethers)
