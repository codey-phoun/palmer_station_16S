#!/bin/bash
file=16S_libraries_combined_taxonomy_feature_table_unassigned
blastn -db ~/blastdb/16S_ribosomal_RNA \
-query ${file}.fasta \
-outfmt "10 qseqid sseqid length qcovs pident evalue staxids stitle" \
-task blastn \
-max_target_seqs 10 \
> ./blast_results/${file}_blast_results.csv

# add header to csv
sed -i 1i"qseqid,sseqid,length,qcovs,pident,evalue,staxids,stitle" ./blast_results/${file}_blast_results.csv