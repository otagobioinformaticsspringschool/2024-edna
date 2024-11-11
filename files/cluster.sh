vsearch --cluster_size sorted_combined.fasta --centroids centroids.fasta --sizein --id 0.97 --sizeout

vsearch --uchime3_denovo centroids.fasta --sizein --fasta_width 0 --nonchimeras otus.fasta --relabel OTU.

vsearch --usearch_global combined.fasta --db otus.fasta --id 0.9 --otutabout otu_frequency_table.tsv