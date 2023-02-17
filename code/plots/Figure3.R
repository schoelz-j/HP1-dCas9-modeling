rm(list = ls())
setwd('/Users/jack/Box/results/revisions/expression/')
read.csv('/Users/jack/Box/results/schoelz_feng_riddle_2020/2019-09-11/DESeq All/YW vs HP1a/KO vs YW DESeq2 annotated results with normalized counts.csv')
hp1a <- read.csv('/Users/jack/Box/results/schoelz_feng_riddle_2020/2019-09-11/DESeq All/YW vs HP1a/KO vs YW DESeq2 annotated results with normalized counts.csv')
hp1b <- read.csv('/Users/jack/Box/results/schoelz_feng_riddle_2020/2019-09-11/DESeq All/YW vs HP1B KO2/KO vs YW DESeq2 annotated results with normalized counts.csv')
hp1c <- read.csv('/Users/jack/Box/results/schoelz_feng_riddle_2020/2019-09-11/DESeq All/YW vs HP1C/KO vs YW DESeq2 annotated results with normalized counts.csv')

library(ggplot2)
library(dplyr)
install.packages('ggrastr')
library(ggrastr)

hp1a <- hp1a %>%
  mutate(Neg.Log.Padj = -log(padj)) %>%
  mutate(significant = case_when(padj < 0.05 ~ 1,
                                 padj > 0.05 ~ 0))
hp1a$significant <- as.factor(hp1a$significant)

hp1b <- hp1b %>%
  mutate(Neg.Log.Padj = -log(padj)) %>%
  mutate(significant = case_when(padj < 0.05 ~ 1,
                                 padj > 0.05 ~ 0))
hp1b$significant <- as.factor(hp1b$significant)

hp1c <- hp1c %>%
  mutate(Neg.Log.Padj = -log(padj)) %>%
  mutate(significant = case_when(padj < 0.05 ~ 1,
                                 padj > 0.05 ~ 0))
hp1c$significant <- as.factor(hp1c$significant)
str(hp1a)
colnames(hp1a)[1] <- 'FBgn'
colnames(hp1b)[1] <- 'FBgn'
colnames(hp1c)[1] <- 'FBgn'

plus_data <- read.csv('/Users/jack/Box/results/modeling_study/datasets/model_inputs/hp1_plus_input.csv')
hp1a_plus <- left_join(plus_data, hp1a, by = 'FBgn')
hp1b_plus <- left_join(plus_data, hp1b, by = 'FBgn')
hp1c_plus <- left_join(plus_data, hp1c, by = 'FBgn')

hp1a_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(GAGA_Sites >= 1), y = log2FoldChange, fill = as.factor(GAGA_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#473E92'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
hp1a_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(DRE_Sites >= 1), y = log2FoldChange, fill = as.factor(DRE_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#473E92'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
hp1a_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(Disco_Sites >= 1), y = log2FoldChange, fill = as.factor(Disco_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#473E92'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')

bgaga <- hp1b_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(GAGA_Sites >= 1), y = log2FoldChange, fill = as.factor(GAGA_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#2D8775'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
bdre <- hp1b_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(DRE_Sites >= 1), y = log2FoldChange, fill = as.factor(DRE_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#2D8775'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
bdisco <- hp1b_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(Disco_Sites >= 1), y = log2FoldChange, fill = as.factor(Disco_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#2D8775'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')

cgaga <- hp1c_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(GAGA_Sites >= 1), y = log2FoldChange, fill = as.factor(GAGA_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#8F3086'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
cdre <- hp1c_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(DRE_Sites >= 1), y = log2FoldChange, fill = as.factor(DRE_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#8F3086'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
cdisco <- hp1c_plus %>%
  filter(significant == 1) %>%
  ggplot(aes(x = as.factor(Disco_Sites >= 1), y = log2FoldChange, fill = as.factor(Disco_Sites >= 1)))+
  geom_boxplot(width = 0.5, notch = T)+
  scale_fill_manual(values = c('grey', '#8F3086'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = 'NA', color = 'black'),
        legend.position = 'none')
hp1b_plus <- hp1b_plus %>%
  mutate(DRE = as.factor(DRE_Sites >= 1))

mean(hp1b_plus$YW1[hp1b_plus$DRE == T & hp1b_plus$significant == 1], na.rm = T)
mean(hp1b_plus$YW2[hp1b_plus$DRE == T & hp1b_plus$significant == 1], na.rm = T)
mean(hp1b_plus$bKO2a[hp1b_plus$DRE == T & hp1b_plus$significant == 1], na.rm = T)
mean(hp1b_plus$bKO2b[hp1b_plus$DRE == T & hp1b_plus$significant == 1], na.rm = T)

pdf('Figure6_d3.pdf', height = 8, width = 9)
ggarrange(bgaga, bdre, bdisco, cgaga, cdre, cdisco,
          ncol = 3, nrow = 2, labels = c('A.', 'B.', 'C.', 'D.', 'E.', 'F.'))
dev.off()
getwd()
nrow(plus_data)

length(hp1b_plus$FBgn[hp1b_plus$significant == 1])
length(hp1b_plus$FBgn[hp1b_plus$DRE_Sites >= 1])
length(hp1b_plus$FBgn[hp1b_plus$significant == 1 & hp1b_plus$DRE_Sites >= 1])
c <- nrow(hp1b_plus)
phyper(213, 1835, 5393-1835, 837, lower.tail = T)
length(hp1c_plus$FBgn[hp1c_plus$significant == 1])
length(hp1c_plus$FBgn[hp1c_plus$DRE_Sites >= 1])
length(hp1c_plus$FBgn[hp1c_plus$significant == 1 & hp1b_plus$DRE_Sites >= 1])
nrow(hp1c_plus)
phyper(412, 1835, 5393-1835, 1207, lower.tail = T)

wilcox.test(hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$GAGA_Sites == 0],
            hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$GAGA_Sites >= 1])
wilcox.test(hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$Disco_Sites == 0],
            hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$Disco_Sites >= 1])
wilcox.test(hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$DRE_Sites == 0],
            hp1b_plus$log2FoldChange[hp1b_plus$significant == 1 & hp1b_plus$DRE_Sites >= 1])

wilcox.test(hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$GAGA_Sites == 0],
            hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$GAGA_Sites >= 1])
wilcox.test(hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$Disco_Sites == 0],
            hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$Disco_Sites >= 1])
wilcox.test(hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$DRE_Sites == 0],
            hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$DRE_Sites >= 1])

str(hp1c)
median(hp1c_plus$log2FoldChange[hp1c_plus$significant == 1 & hp1c_plus$DRE_Sites >= 1], na.rm = T)
mean(hp1c_plus$YW1[hp1c_plus$significant == 1 & hp1c_plus$DRE_Sites >= 1], na.rm = T)
