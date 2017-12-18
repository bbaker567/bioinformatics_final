#This file is for testing. See liz_functions.R for completed functions.
setwd("~/Desktop/Bioinformatics_Final/")


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
  new_df_results$Top_Classification <- NA
  new_df_results$Super_Kingdom <- NA
  new_df_results$Phylum <- NA
  new_df_results$Class <- NA
  new_df_results$Order <- NA
  new_df_results$Family <- NA
  new_df_results$Genus <- NA
  
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

superkingdom_stacked_plot <- function(metagenome_df, super_kingdom, division, main_title, y_label){
  library(ggplot2)
  
  only_super_kingdom <- subset(metagenome_df, Super_Kingdom == super_kingdom)
  abundancy_df <- as.data.frame(sort(table(only_super_kingdom[, division]), decreasing = F))
  abundancy_df$Freq = 0
  
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

superkingdom_pie_chart(QQ_1_to_10, "Microbe Respresentation in QQ Hot Spring by Superkingdom")

condense_species <- function(metagenome_df){
  ordered_metagenome_df <- metagenome_df[order(metagenome_df$Species), ]
  ordered_metagenome_df <- species_name_cleanup(ordered_metagenome_df)
  df_loop_length <- nrow(ordered_metagenome_df)
  
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

test_plot <- superkingdom_stacked_plot(QQ_1_to_10, "Bacteria", "Phylum", "Representation of Bacterial Phyla in QQ Hot Spring", "Phyla Abundance")
test_plot

QQ_species_condensed <- condense_species(new_QQ_1_to_10)

?factor
bacteria_phyla_stacked_plot(QQ_1_to_10, "Representation of Bacterial Phyla in QQ Hot Spring", "Phyla Abundance")
QQ_1_to_10 <- metagenome_sample_blast("QQ_1_to_10.csv", 301, .7)
result_sum <- sum(QQ_1_to_10$Frequency)
result_sum

cells_only <- subset(QQ_1_to_9, Top_Classification == "cellular organisms")
pie(table(cells_only$Super_Kingdom))

bacteria <- subset(QQ_1_to_9, Super_Kingdom == "Bacteria")
pie(table(bacteria$Phylum))

archaea <- subset(QQ_1_to_9, Super_Kingdom == "Archaea")
pie(table(archaea$Phylum))

bacteria_phyla <- as.data.frame(table(bacteria$Phylum))
bacteria_ordered <- bacteria_phyla[order(bacteria_phyla$Freq), ]
stacked_bar_plot <- ggplot(bacteria_phyla[order(bacteria_phyla$Freq), ], aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")

stacked_bar_plot 

bacteria_phyla$Var1 <- factor(bacteria_phyla$Var1, levels = rev(levels(bacteria_phyla$Var1)))
ggplot(bacteria_phyla, aes(fill = Var1, order = -as.numeric(Freq))) 
  + geom_bar()

ggplot(bacteria_phyla, aes(Var1, Freq, fill=Var1)) + geom_bar(stat="identity")

ggplot(bacteria_phyla, aes(Var1, Freq, group = Freq, stat='histogram')) +
  geom_(aes(fill = Var1)) +
  coord_flip()

ggplot((bacteria_phyla), aes(x=Var1, y=Freq, fill=Var1))+
  geom_bar(stat="identity") + labs(title="Original dataframe")

ggplot(bacteria_ordered[order(bacteria_phyla$Freq), ], aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity") + labs(title="Original dataframe")

ggplot(test, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity")

test <- bacteria_ordered[order(-bacteria_phyla$Freq),]
test

test_2 <- as.data.frame(sort(table(bacteria$Phylum), decreasing = F))

ggplot(test_2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity") + labs(title="Original dataframe")

x=1
new_QQ_1_to_10[x, "Species"] == new_QQ_1_to_10[x + 1, "Species"]
new_QQ_1_to_10[-(x + 1), ]

new_QQ_1_to_10 <- QQ_1_to_10
new_QQ_1_to_10 <- QQ_1_to_10[order(QQ_1_to_10$Species), ]
new_QQ_1_to_10 <- species_name_cleanup(new_QQ_1_to_10)
df_loop_length <- nrow(new_QQ_1_to_10)
test <- condense_species(new_QQ_1_to_10)

for(x in 1:df_loop_length){
  if(is.na(new_QQ_1_to_10[x, "Species"] == new_QQ_1_to_10[(x + 1), "Species"])){
    break
  }
  while(new_QQ_1_to_10[x, "Species"] == new_QQ_1_to_10[(x + 1), "Species"]){
    new_QQ_1_to_10[x, "Frequency"] = new_QQ_1_to_10[x, "Frequency"] + new_QQ_1_to_10[x + 1, "Frequency"]
    new_QQ_1_to_10 <- new_QQ_1_to_10[-(x + 1), ]
  }
}

new_QQ_1_to_10[30, "Species"] == new_QQ_1_to_10[33, "Species"]

summed <- rowSums(zscore[, c(1, 2, 3, 5)])


bacteria_phyla_stacked_plot <- function(metagenome_df, main_title, y_label){
  library(ggplot2)
  
  bacteria <- subset(metagenome_df, Super_Kingdom == "Bacteria")
  sorted_table <- as.data.frame(sort(table(bacteria$Phylum), decreasing = F))
  
  stacked_bar_plot <- ggplot(sorted_table, aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity") + labs(title=main_title) +
    xlab("") + ylab(y_label)
  
  return(stacked_bar_plot)
}

bacteria <- subset(new_QQ_1_to_10, Super_Kingdom == "Bacteria")
sorted_table <- as.data.frame(sort(table(bacteria[, "Phylum"]), decreasing = F))

division_names <- sorted_table$Var1
abundancy_df <- as.data.frame(division_names)
names(abundancy_df)[names(abundancy_df) == '`sorted_table$Var1`'] <- 'Division'
abundancy_df$Count <- 0

for(x in 1:nrow(bacteria)){
  phylum <- bacteria[x, "Phylum"]
  abundancy_row <- which(abundancy_df$division_names == phylum)
  abundancy_df[abundancy_row, "Count"] <- abundancy_df[abundancy_row, "Count"] + bacteria[x, "Frequency"]
}

stacked_bar_plot <- ggplot(sorted_table, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") + labs(title="main_title") +
  xlab("") + ylab("y_label")
