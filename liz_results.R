setwd("~/Desktop/Bioinformatics_Final/")
library(ape)
library(seqinr)

genome_results <- read.csv(file = "results_QQ_subset_500_R_2.csv", header = F)
colnames(genome_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

high_coverage <- subset(genome_results, length >= 210)
tabled_results <- table(high_coverage$sseqid)
df_results <- as.data.frame(tabled_results, stringsAsFactors=FALSE)
colnames(df_results) <- c("Accession", "Frequency")
df_results$Species <- NA
new_df_results <- df_results[-c(which(df_results$Frequency == 0)), ] 
result_length <- nrow(new_df_results)

for(x in 1:result_length){
  accession <- new_df_results[x, "Accession"]
  sequence <- read.GenBank(accession)
  species <- attr(sequence, "species")
  new_df_results[x, "Species"] <- species
}

sum(new_df_results$Frequency, na.rm = T)

?subset
