fasta_subsample <- function(big_fasta, new_length, new_file){
  library(seqinr)
  long_reads <- read.fasta(big_fasta, as.string = T)
  subset <- sample(long_reads, new_length)
  subset_unlist <- unlist(subset)
  
  write.fasta(sequences = subset, names = names(subset_unlist), file.out = new_file)
}

metagenome_sample_blast <- function(csv_file, read_length, minimum_coverage){
  library(ape)
  library(seqinr)
  library(taxize)
  
  genome_results <- read.csv(csv_file, header = F)
  colnames(genome_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  high_coverage <- subset(genome_results, length >= (minimum_coverage*read_length))
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
    species <- gsub("-", " ", species)
    species <- gsub("\\[", "", species)
    species <- gsub("]", "", species)
    new_df_results[x, "Species"] <- species
    
    domain_list <- classification(species, db = 'ncbi')
    unlisted_domain <- unlist(domain_list)
    new_df_results[x, "Domain"] <- unlisted_domain[2]
  }
  
  return(new_df_results)
}

QQ_1_to_6 <- metagenome_sample_blast("subsets_1_to_6.csv", 301, .7)

pie(new_df_results$Frequency, labels = new_df_results$Species)

pie(table(new_df_results$Domain))

samples_blasted <- 500
QQ_read_total <- 4000000
total_sample_proportion <- 500/4000000
total_sample_identified <- sum(new_df_results$Frequency, na.rm = T)
id_total_sample_proportion <- (total_sample_identified / samples_blasted) * 100