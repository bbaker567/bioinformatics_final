#Read in blast results file, add column headers
QQ_DF_raw <- read.csv("QQ_blastn_results_brit", sep = "\t", header = F, 
                    col.names = c("qseqid", "sseqid", "pident", "length", 
                                  "mismatch", "gapopen", "qstart", "qend", "sstart", 
                                  "send", "evaule", "bitscore", "qcovs"))

#Only keep results with a % length of >= 70%
QQ_DF_subest <- subset(QQ_DF_raw, QQ_DF_raw$length >=210)
sseqid <- QQ_DF_subest$sseqid

#Plot of number of hits per species
library(ggplot2)
ggplot(data = QQ_DF_subest) +
  geom_bar(mapping = aes(x = sseqid, fill = sseqid), colour = "black") +
  coord_flip() +
 theme(legend.position = "bottom") +
  ylab("Number of blast Hits") +
  xlab("Species") + 
  labs(title = "Number of blast hits per Species") 
