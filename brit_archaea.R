#Read in .csv file, add headers
QQ_DF_raw <- read.csv("QQ_blastn_results_brit", sep = "\t", header = F, 
                    col.names = c("qseqid", "sseqid", "pident", "length", 
                                  "mismatch", "gapopen", "qstart", "qend", "sstart", 
                                  "send", "evaule", "bitscore", "qcovs"))

#Only keep results with a % length of >= 70%
QQ_DF_subest <- subset(QQ_DF_raw, QQ_DF_raw$length >=210)
sseqid <- QQ_DF_subest$sseqid

#Makes a data frame of the frequency of Genus hits
tabled_results <- table( QQ_DF_subest$sseqid)
genus_frequency <- as.data.frame(tabled_results, stringsAsFactors=FALSE)

#Plot of number of hits per species
library(ggplot2)
ggplot(data = QQ_DF_subest) +
  geom_bar(mapping = aes(x = sseqid, fill = sseqid), colour = "black") +
  coord_flip() +
  ylab("Number of blast Hits") +
  xlab("Genus") + 
  labs(title = "Number of BLAST hits per Archaeal Genus") +
  theme(legend.position="none")

#Extract the frequency of just the species by e-value
install.packages("ape")
library(ape)
QQ_DF_sseqid <- ddply(QQ_DF_1, .(sseqid), c("nrow"))

#Look at the ncbi pyhlum classification
install.packages("taxize")
library(taxize)
test <- classification(s_1, db= "ncbi")
test_1 <- tax_name(s_1, get = "phylum", db = "ncbi")

#Make a datframe based off of taxize results 
df_subset <- QQ_DF_1[1:5,]
s_1 <- df_subset$sseqid 
QQ_DF_raw <- read.csv("QQ_blastn_results_brit", sep = "\t", header = F, 
                      col.names = c("qseqid", "sseqid", "pident", "length", 
                                    "mismatch", "gapopen", "qstart", "qend", "sstart", 
                                    "send", "evaule", "bitscore", "qcovs"))

QQ_DF_subest <- subset(QQ_DF_raw, QQ_DF_raw$length >=210)
sseqid <- QQ_DF_subest$sseqid
tabled_results <- table(QQ_DF_subest$sseqid)
df_results <- as.data.frame(tabled_results, stringsAsFactors=FALSE)
colnames(df_results) <-  c("Genus", "Frequency")
new_df_results <- df_results[-c(which(df_results$Frequency == 0)), ] 
new_df_results$Genus -> Genus
full_classification <- classification(Genus, db = 'ncbi')
unlisted_classification <- unlist(full_classification)
new_df_results[, "Top_Classification"] <- unlisted_classification[1]
new_df_results[, "Super_Kingdom"] <- unlisted_classification[2]
new_df_results[, "Phylum"] <- unlisted_classification[3]
new_df_results[, "Class"] <- unlisted_classification[4]
new_df_results[, "Order"] <- unlisted_classification[5]
new_df_results[, "Family"] <- unlisted_classification[6]
new_df_results[, "Genus"] <- unlisted_classification[7]

#Make a plot based on the phylum
ggplot(data = Phylum) +
  geom_bar(mapping = aes(x = Phylum_1), colour = "black") +
  coord_flip() +
  theme(legend.position = "bottom") +
  ylab("Number of blast Hits") +
  xlab("Species") + 
  labs(title = "Number of blast hits per Species") 
