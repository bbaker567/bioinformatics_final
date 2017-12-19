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
