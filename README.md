# bioinformatics_final
Final Project

#Step 1
##Collected genomes from NCBI as refernce genomes that we thought could potentially be in the QQ metagenome, and made a blast #database.

makeblastdb -dbtype nucl -in bunchofeditedgenomes.fna -out dbbunchofeditedgenomes

#Step 2
##Converted our fastq file to a fasta file 

sed -n '1~4s/^@/>/p;2~4p' QQ_small.fastq > QQ_small.fasta

#Step 3
##Ran a blastn search of the reference genomes against the QQ metagenome fasta file. We limitied the results to only the top #hit, and set the e value to 1e-5

blastn -query QQ_metagenome.fasta -db dbbunchofeditedgenomes -evalue 1e-5 -outfmt 6 -max_target_seqs 1 >> new_results
