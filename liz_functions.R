#this reads a FASTA file, takes a random subset of its sequences, and writes it to a new FASTA file
#enter the file name, the desired number of sequences to sample, and a new file name
fasta_subsample <- function(big_fasta, new_length, new_file){
  library(seqinr)
  long_reads <- read.fasta(big_fasta, as.string = T)
  
  #This takes a random sample of the reads in the file, according to initial parameters
  subset <- sample(long_reads, new_length)
  subset_unlist <- unlist(subset)
  
  #This outputs a new FASTA file with the sampled reads
  write.fasta(sequences = subset, names = names(subset_unlist), file.out = new_file)
}

#This function reads a CSV file with columns arising from standard tabular output (outfmt 6) 
#from command-line blastn.
#It returns a data frame with columns for accession number, frequency,
#species, and further columns for taxonomic groupings.
metagenome_sample_blast <- function(csv_file, read_length, minimum_coverage){
  library(ape)
  library(seqinr)
  library(taxize)
  
  genome_results <- read.csv(csv_file, header = F)
  colnames(genome_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  #This eliminates matches with low coverage, which are of lower quality than long coverage matches.
  high_coverage <- subset(genome_results, length >= (minimum_coverage*read_length))
  tabled_results <- table(high_coverage$sseqid)
  
  df_results <- as.data.frame(tabled_results, stringsAsFactors=FALSE)
  colnames(df_results) <- c("Accession", "Frequency")
  df_results$Species <- NA
  
  new_df_results <- df_results[-c(which(df_results$Frequency == 0)), ] 
  result_length <- nrow(new_df_results)
  new_df_results$Top_Classification <- NA
  new_df_results$Super_Kingdom <- NA
  new_df_results$Phylum <- NA
  new_df_results$Class <- NA
  new_df_results$Order <- NA
  new_df_results$Family <- NA
  new_df_results$Genus <- NA
  
  #This acquires the match's accession number, looks up the species name, and stores 
  #phylogenetic information about each species.
  for(x in 1:result_length){
    accession <- new_df_results[x, "Accession"]
    sequence <- read.GenBank(accession)
    species <- attr(sequence, "species")
    species <- gsub("\\[", "", species)
    species <- gsub("]", "", species)
    new_df_results[x, "Species"] <- species
    
    full_classification <- classification(species, db = 'ncbi')
    unlisted_classification <- unlist(full_classification)
    new_df_results[x, "Top_Classification"] <- unlisted_classification[1]
    new_df_results[x, "Super_Kingdom"] <- unlisted_classification[2]
    new_df_results[x, "Phylum"] <- unlisted_classification[3]
    new_df_results[x, "Class"] <- unlisted_classification[4]
    new_df_results[x, "Order"] <- unlisted_classification[5]
    new_df_results[x, "Family"] <- unlisted_classification[6]
    new_df_results[x, "Genus"] <- unlisted_classification[7]
  }
  
  return(new_df_results)
}

#This function cleans up and simplifies the names of species.
species_name_cleanup <- function(formatted_df){
  df_length <- nrow(formatted_df)
  
  for(x in 1:df_length){
    species <- formatted_df[x, "Species"]
    species <- gsub("_", " ", species)
    #species <- gsub("\\s[[:graph:]]*[0-9]+.*", "", species)
    #species <- gsub("\\ssp\\.", "", species)
    #species <- gsub("\\s[A-Z]{2,}", "", species)
    
    formatted_df[x, "Species"] <- species
  }
  
  return(formatted_df)
}

#This cleans up the dataframe by combining rows with the same species name that may have arisen 
#due to BLAST returning sequences with different accession numbers from the same species.
condense_species <- function(metagenome_df){
  ordered_metagenome_df <- metagenome_df[order(metagenome_df$Species), ]
  ordered_metagenome_df <- species_name_cleanup(ordered_metagenome_df)
  df_loop_length <- nrow(ordered_metagenome_df)
  
  #This preserves the frequency counts of species, regardless of how many different accession numbers 
  #showed up per species. 
  for(x in 1:df_loop_length){
    if(is.na(ordered_metagenome_df[x, "Species"] == ordered_metagenome_df[(x + 1), "Species"])){
      break
    }
    while(ordered_metagenome_df[x, "Species"] == ordered_metagenome_df[(x + 1), "Species"]){
      ordered_metagenome_df[x, "Frequency"] = ordered_metagenome_df[x, "Frequency"] + ordered_metagenome_df[x + 1, "Frequency"]
      ordered_metagenome_df <- ordered_metagenome_df[-(x + 1), ]
    }
  }
  return(ordered_metagenome_df)
}

#This function creates a pie chart showing the number of matches per super kingdom.
#Pie charts are as useful for more varied data, but the differences in super kingdom variation
#are so stark that a pie chart paints a clear picture.
superkingdom_pie_chart <- function(metagenome_df, main_title){
  cells_only <- subset(metagenome_df, Top_Classification == "cellular organisms")
  
  abundancy_df <- as.data.frame(sort(table(cells_only[, "Super_Kingdom"]), decreasing = F))
  abundancy_df$Freq = 0
  
  for(x in 1:nrow(cells_only)){
    one_division <- cells_only[x, "Super_Kingdom"]
    abundancy_row <- which(abundancy_df$Var1 == one_division)
    abundancy_df[abundancy_row, "Freq"] <- abundancy_df[abundancy_row, "Freq"] + cells_only[x, "Frequency"]
  }
  
  abundancy_df$Var1 <- factor(abundancy_df$Var1, levels = abundancy_df$Var1[order(rev(abundancy_df$Freq))])
  
  stacked_bar_plot <- ggplot(abundancy_df, aes(x="", y=Freq, fill=Var1, label = Freq)) +
    geom_bar(width = 1, stat = "identity") + labs(title=main_title) +
    xlab("") + ylab("")+
    coord_polar(theta = "y") +
    geom_text(size = 3, position = position_stack(vjust = 0.5))
  
  return(stacked_bar_plot)
}

#This function produces a stacked barplot for a given super kingdom and division, such as Phylum or Order.
#The height of the bars refers to the number of total hits per given division, not the number of species
#within it.
superkingdom_stacked_plot <- function(metagenome_df, super_kingdom, division, main_title, y_label){
  library(ggplot2)
  
  only_super_kingdom <- subset(metagenome_df, Super_Kingdom == super_kingdom)
  abundancy_df <- as.data.frame(sort(table(only_super_kingdom[, division]), decreasing = F))
  
  for(x in 1:nrow(only_super_kingdom)){
    one_division <- only_super_kingdom[x, division]
    abundancy_row <- which(abundancy_df$Var1 == one_division)
    abundancy_df[abundancy_row, "Freq"] <- abundancy_df[abundancy_row, "Freq"] + only_super_kingdom[x, "Frequency"]
  }
  
  abundancy_df$Var1 <- factor(abundancy_df$Var1, levels = abundancy_df$Var1[order(abundancy_df$Freq)])
  
  stacked_bar_plot <- ggplot(abundancy_df, aes(x="", y=Freq, fill=Var1, order=Var1)) +
    geom_bar(width = 1, stat = "identity") + labs(title=main_title) +
    xlab("") + ylab(y_label)
  
  return(stacked_bar_plot)
}