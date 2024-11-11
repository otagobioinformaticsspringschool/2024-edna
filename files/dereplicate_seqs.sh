source ~/obss_2021/edna/scripts/envs.sh

vsearch --derep_fulllength ../data/fasta/combined.fasta --minuniquesize 2 --sizeout --output derep_combined.fasta --relabel Uniq.

vsearch --sortbysize derep_combined.fasta --output sorted_combined.fasta