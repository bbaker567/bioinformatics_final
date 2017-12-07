fasta_subsample <- function(big_fasta, new_length, new_file){
  library(seqinr)
  long_reads <- read.fasta(big_fasta, as.string = T)
  subset <- sample(long_reads, new_length)
  subset_unlist <- unlist(subset)
  
  write.fasta(sequences = subset, names = names(subset_unlist), file.out = new_file)
}