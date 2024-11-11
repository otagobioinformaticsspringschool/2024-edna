module load eDNA

cd ~/obss_2021/edna/

cutadaptQC \
 -m ./docs/sample_metadata.tsv \
 -r ./data/FTP103_S1_L001_R1_001.fastq.gz \
 -f ./data \
 -t 4 \
 -n fish_project