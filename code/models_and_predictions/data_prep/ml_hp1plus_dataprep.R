# Integrate data for HP1+genomic features ml model
setwd('/Users/jack/Box/results/gene_density/')

library(dplyr)
library(ggplot2)
# need one_hot function
#install.packages('mltools')
library(mltools)
library(data.table)

# HP1 Coverage at TSS data
# (Also includes target variable)
hp1 <- read.csv('data/HP1-coverage/s2/HP1_ml_data_8_26_2021.csv')
str(hp1)
# 9528 genes

# Previous dataset - GC Body + Promoter, Gene density, strand,
# Pausing index
# Remove non-finite values of PInd
prev <- read.csv('model_data_8_19_2021.csv')
str(prev)
prev <- prev %>%
  select(Density, GC.Promoter, GC.Body, Strand, PInd, FBgn) %>%
  filter(is.finite(PInd) == T)
# 11739

# Chromatin State
cstate <- read.csv('S2_Chromatin_States_RF_Input.csv')
str(cstate)
cstate$state.name <- as.factor(cstate$state.name)
cstate_1hot <- one_hot(as.data.table(cstate))
# 17253

# GAGA motif
gaga_all <- read.csv('gaga_fimo_9_2_2021.tsv', sep = '\t')
str(gaga_all)

colnames(gaga_all)[3] <- 'FBgn'
gaga <- gaga_all %>%
  group_by(FBgn) %>%
  summarise(GAGA_Sites = n(),
            GAGA_Avg_Score = mean(score))

# DRE motif
dre_all <- read.csv('dre_fimo_results.tsv', sep = '\t')
str(dre_all)

colnames(dre_all)[3] <- 'FBgn'
dre <- dre_all %>%
  group_by(FBgn) %>%
  summarise(DRE_Sites = n(),
            DRE_Avg_Score = mean(score))

# Disco motif
disco_all <- read.csv('disco_fimo_results.tsv', sep = '\t')

colnames(disco_all)[3] <- 'FBgn'
disco <- disco_all %>%
  group_by(FBgn) %>%
  summarise(Disco_Sites = n(),
            Disco_Avg_Score = mean(score))

# Exon number
exons <- read.csv('drosophila_exon_numbers.csv')
str(exons)

# Nucleotide Composition
nucs <- read.csv('promoter_nucleotide_frequencies.csv')
str(nucs)

# TSS accessibility
tss_acc <- read.csv('data/atac/TSS_Accessibility_S2.csv')
str(tss_acc)

# Integrate datasets
ml_plus <- inner_join(hp1, prev, by = 'FBgn')
ml_plus <- inner_join(ml_plus, cstate_1hot, by = 'FBgn')
ml_plus <- inner_join(ml_plus, exons, by = 'FBgn')
ml_plus <- inner_join(ml_plus, nucs, by = 'FBgn')
ml_plus <- inner_join(ml_plus, tss_acc, by = 'FBgn')

# For motifs, want values of 0 for genes with missing data
ml_plus <- left_join(ml_plus, dre)
ml_plus <- left_join(ml_plus, gaga)
ml_plus <- left_join(ml_plus, disco)

# Set NA values in these motifs to 0
ml_plus$DRE_Sites[is.na(ml_plus$DRE_Sites)] <- 0
ml_plus$DRE_Avg_Score[is.na(ml_plus$DRE_Avg_Score)] <- 0
ml_plus$Disco_Sites[is.na(ml_plus$Disco_Sites)] <- 0
ml_plus$Disco_Avg_Score[is.na(ml_plus$Disco_Avg_Score)] <- 0
ml_plus$GAGA_Sites[is.na(ml_plus$GAGA_Sites)] <- 0
ml_plus$GAGA_Avg_Score[is.na(ml_plus$GAGA_Avg_Score)] <- 0



#

# Output data for scikit learn modeling

write.table(ml_plus, file = 'hp1-plus-model-data.csv',
            col.names = T, row.names = F, quote = F, sep = ',')


no_hp1 <- ml_plus %>%
  select(-HP1a, -HP1B, -HP1C, -A_B, -A_C, -B_C, -A_B_C)
write.table(no_hp1, file = 'no-hp1-model-data.csv', col.names = T,
            row.names = F, quote = F, sep = ',')

# Generate 'permuted' random datasets to test robustness
library(picante)

ml_mat <- ml_plus %>%
  select(-FBgn, -log_S2) %>%
  as.matrix()
mlplus_random <- randomizeMatrix(as.matrix(ml_mat), null.model = 'richness',
                iterations = 10000)
no_mat <- no_hp1 %>%
  select(-FBgn, -log_S2) %>%
  as.matrix()
no_random <- randomizeMatrix(as.matrix(no_mat), null.model = 'richness',
                             iterations = 10000)

hp1_mat <- ml_plus %>%
  select(HP1a, HP1B, HP1C, A_B, A_C, B_C, A_B_C, -FBgn, -log_S2) %>%
  as.matrix()
hp1_random <- randomizeMatrix(as.matrix(hp1_mat), null.model = 'richness',
                              iterations = 10000)

mlplus_random <- data.frame(mlplus_random)
mlplus_random$log_S2 <- ml_plus$log_S2
no_random <- data.frame(no_random)
no_random$log_S2 <- no_hp1$log_S2
hp1_random <- data.frame(hp1_random)
hp1_random$log_S2 <- ml_plus$log_S2


write.table(mlplus_random, file = 'HP1plus_randomized.csv', col.names = T, row.names = F, sep = ',',
            quote = F)
write.table(no_random, file = 'NoHP1_randomized.csv', col.names = T, row.names = F, sep = ',',
            quote = F)
write.table(hp1_random, file = 'HP1only_randomized.csv', col.names = T, row.names = F, sep = ',',
            quote = F)
