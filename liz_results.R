setwd("~/Desktop/Bioinformatics_Final/")
library(ape)
library(seqinr)
library(taxize)

genome_results <- read.csv(file = "subsets_1_to_6.csv", header = F)
colnames(genome_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

high_coverage <- subset(genome_results, length >= 210)
tabled_results <- table(high_coverage$sseqid)
df_results <- as.data.frame(tabled_results, stringsAsFactors=FALSE)
colnames(df_results) <- c("Accession", "Frequency")
df_results$Species <- NA
new_df_results <- df_results[-c(which(df_results$Frequency == 0)), ] 
result_length <- nrow(new_df_results)
new_df_results$Domain <- NA

for(x in 1:result_length){
  accession <- new_df_results[x, "Accession"]
  sequence <- read.GenBank(accession)
  species <- attr(sequence, "species")
  species <- gsub("_", " ", species)
  species <- gsub("\\[", "", species)
  species <- gsub("]", "", species)
  new_df_results[x, "Species"] <- species
  
  domain_list <- classification(species, db = 'ncbi')
  unlisted_domain <- unlist(domain_list)
  new_df_results[x, "Domain"] <- unlisted_domain[2]
}

pie(new_df_results$Frequency, labels = new_df_results$Species)

pie(table(new_df_results$Domain))

samples_blasted <- 500
QQ_read_total <- 4000000
total_sample_proportion <- 500/4000000
total_sample_identified <- sum(new_df_results$Frequency, na.rm = T)
id_total_sample_proportion <- (total_sample_identified / samples_blasted) * 100