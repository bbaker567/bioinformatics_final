setwd("~/Desktop/Bioinformatics_Final/")
fasta_subsample("QQ_metagenome.fna", 1000, "QQ_subset_1000_R_2.fna")
genome_results <- read.csv(file = "59_db_results.tsv", header = T)

high_coverage <- subset(x= genome_results, genome_results$qcovs >= 70)
tabled_results <- table(high_coverage$sseqid)
df_results <- as.data.frame(tabled_results, stringsAsFactors=FALSE)
df_results$Classification <- NA
colnames(df_results)[1] <- "Organism"
colnames(df_results)[2] <- "Count"
assigned_sum <- sum(df_results$Count)
df_results[nrow(df_results) + 1, 1] <- "Unassigned"
df_results[60, 2] <- 4000000 - assigned_sum
pie(df_results$Count, labels = df_results$Organism)
