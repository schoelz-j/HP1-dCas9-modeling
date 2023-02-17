# Calculate number of exons per gene

setwd('/Users/johnschoelz/Box/results/gene_density/')

library(ggplot2)
library(dplyr)

exons <- read.csv('exon.gff', sep = ' ', header = F)
str(exons)
colnames(exons)[2] <- 'FBgn_'

# This is all unnecessary - each line in the file is an exon
#get_fbgn <- function(i){
#  fbgn_ <- strsplit(exons$meta, split = ' ')[[i]][2]
#  substr(fbgn_, 1, nchar(fbgn_)-1)
#}
#exons$fbgn <- unlist(lapply(seq(1:nrow(exons)), get_fbgn))

exons <- exons %>%
  mutate(FBgn = substr(FBgn_, 1, nchar(FBgn_)-1))
exon_num <- exons %>%
  group_by(FBgn) %>%
  summarise(exon_num = n())

write.table(exon_num, file = 'drosophila_exon_numbers.csv', col.names = T,
            row.names = F, quote = F, sep = ',')
