setwd("~/Desktop/Bioinformatics_Final/")
fasta_subsample("QQ_metagenome.fna", 1000, "QQ_subset_1000_R_2.fna")
genome_results <- read.csv(file = "59_db_results.tsv", header = T)

high_coverage <- subset(x= genome_results, genome_results$qcovs >= 70)
pie(table(high_coverage$sseqid))
table(high_coverage$sseqid)
