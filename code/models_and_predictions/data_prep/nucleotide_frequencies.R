# Calculate nucleotide and dinucleotide frequencies for Drosophila promoters

setwd('/Users/jack/Box/results/gene_density/')

library(BiocManager)
#BiocManager::install('Biostrings')
library(Biostrings)

dmel_promoters <- readDNAStringSet('drosophila_promoters.fa')
dinuc_freqs <- dinucleotideFrequency(dmel_promoters)
nuc_freqs <- oligonucleotideFrequency(dmel_promoters, width = 1)
trinuc_freqs <- trinucleotideFrequency(dmel_promoters)

remove_semicolon <- function(i){
  string_x <- genes[i]
  if(identical(';', substring(string_x, nchar(string_x), nchar(string_x)))){
    new <- substr(string_x, 1, nchar(string_x) - 1)
    return(new)
  }
  string_x
}

genes2 <- unlist(lapply(seq(1:length(genes)), remove_semicolon))

# Write out data with nucleotide frequencies
all_frequencies <- data.frame(cbind(nuc_freqs, dinuc_freqs, trinuc_freqs))
all_frequencies$FBgn <- genes2
str(all_frequencies)

write.table(all_frequencies, file = 'promoter_nucleotide_frequencies.csv', row.names = F, col.names = T,
            quote = F, sep = ',')
