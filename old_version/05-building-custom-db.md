---
title: "Creating A Reference Database"
teaching: 15
exercises: 30
questions:
- "What reference databases are available for metabarcoding"
- "How do I make a custom database?"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

## Reference databases

As we have previously discussed, eDNA metabarcoding is the process whereby we amplify a specific gene region of the taxonomic group of interest from an environmental sample. Once the DNA is amplified, we sequence it and up till now, we have been filtering the data based on quality and assigned sequences to samples. The next step in the bioinformatic process will be to assign a taxonomy to each of the sequences, so that we can determine what species are detected through our eDNA analysis.

There are multiple ways to assign taxonomy to a sequence. While this workshop won't go into too much detail about the different taxonomy assignment methods, any method used will require a reference database.

Several reference databases are available online. The most notable ones being NCBI, EMBL, BOLD, SILVA, and RDP. However, recent research has indicated a need to use custom curated reference databases to increase the accuracy of taxonomy assignment. While certain reference databases are available (RDP: 16S microbial; MIDORI: COI eukaryotes; etc.), this necessary data source is missing in most instances.

Hugh and I, therefore, developed a python software program that can generate custom curated reference databases for you, called CRABS: Creating Reference databases for Amplicon-Based Sequencing. While we will be using it here in the workshop, it is unfortunately not yet publicly available, since we are currently writing up the manuscript. However, we hope that we can release the program by the end of this year. We hope this software program will provide you with the necessary tools to help you generate a custom reference database for your specific study!

## CRABS - a new program to build custom reference databases

The CRABS workflow consists of seven parts, including: 

- (i) downloading and importing data from online repositories into CRABS; 
- (ii) retrieving amplicon regions through *in silico* PCR analysis; 
- (iii) retrieving amplicons with missing primer-binding regions through pairwise global alignments using results from *in silico* PCR analysis as seed sequences; 
- (iv) generating the taxonomic lineage for sequences; 
- (v) curating the database via multiple filtering parameters; 
- (vi) post-processing functions and visualizations to provide a summary overview of the final reference database; and 
- (vii) saving the database in various formats covering most format requirements of taxonomic assignment tools.

Today, we will build a custom curated reference database of fish sequences from the 16S gene region. Once we've completed this reference database, you could try building your own shark reference database from the 16S gene region or try the taxonomy assignment of your OTUs with this newly generated reference database to see how it performs.

Before we get started, we will quickly go over the help information of the program. As before, we will first have to source the `eDNA.sh` script to load the program.

```bash
source eDNA.sh
crabs_v1.0.1 -h
```

```
usage: crabs_v1.0.1 [-h] {db_download,db_import,db_merge,insilico_pcr,pga,assign_tax,dereplicate,seq_cleanup,geo_cleanup,visualization,tax_format} ...

creating a curated reference database

positional arguments:
  {db_download,db_import,db_merge,insilico_pcr,pga,assign_tax,dereplicate,seq_cleanup,geo_cleanup,visualization,tax_format}

optional arguments:
  -h, --help            show this help message and exit
```
{: .output}

As you can see, the help documentation shows you how to use the program, plus displays all the 11 functions that are currently available within CRABS. If you would want to display the help documentation of a function, we can type the following command:

```bash
crabs_v1.0.1 db_download -h
```

```
usage: crabs_v1.0.1 db_download [-h] -s SOURCE [-db DATABASE] [-q QUERY] [-o OUTPUT] [-k ORIG] [-e EMAIL] [-b BATCHSIZE]

downloading sequence data from online databases

optional arguments:
  -h, --help            show this help message and exit
  -s SOURCE, --source SOURCE
                        specify online database used to download sequences. Currently supported options are: (1) ncbi, (2) embl, (3) mitofish, (4) bold, (5)
                        taxonomy
  -db DATABASE, --database DATABASE
                        specific database used to download sequences. Example NCBI: nucleotide. Example EMBL: mam*. Example BOLD: Actinopterygii
  -q QUERY, --query QUERY
                        NCBI query search to limit portion of database to be downloaded. Example: "16S[All Fields] AND ("1"[SLEN] : "50000"[SLEN])"
  -o OUTPUT, --output OUTPUT
                        output file name
  -k ORIG, --keep_original ORIG
                        keep original downloaded file, default = "no"
  -e EMAIL, --email EMAIL
                        email address to connect to NCBI servers
  -b BATCHSIZE, --batchsize BATCHSIZE
                        number of sequences downloaded from NCBI per iteration. Default = 5000
```
{: .output}


> ## Step 1: download sequencing data from online repositories
>During this workshop, you will not have to run this command, since the file has been provided to you!
>
>When we look at the help documentation of the `db_download` function, we can see that CRABS let's us download sequencing data from 4 repositories using the `--source` parameter, including:
>
> - NCBI
> - BOLD
> - EMBL
> - MitoFish
>
>Today, we will be downloading 16S fish sequences from NCBI. According to the help documentation, we will need to fill out multiple parameters.
>
>`--source`: online repository to download sequencing data from.
>
>`--database`: sequencing database of NCBI to download sequences from.
>
>`--query`: NCBI query search that encompasses the data you would like to download. This query is based on the `search details` box on the NCBI website (https://www.ncbi.nlm.nih.gov/nuccore). Let's have a look at this now!
>
>`--output`: output filename.
>
>`--email`: email to let the NCBI servers know who you are. This is a requirement!
>
>The code to download the 16S sequences of sharks would look for example like this:
>
>```bash
>crabs_v1.0.1 db_download \
>--source ncbi \
>--database nucleotide \
>--query '16S[All Fields] AND ("Chondrichthyes"[Organism] OR Chondrichthyes[All Fields]) AND ("1"[SLEN] : "50000"[SLEN])' \
>--output 16S_chondrichthyes_ncbi.fasta \
>--email johndoe@gmail.com
>```
{: .challenge}

Rather than everyone downloading sequencing data off NCBI at the same time, we will copy a sequencing file from the `resources/day` folder in our `references` folder. This file named `ncbi_16S.fasta` is the file generated by the `db_download` function and contains all 16S sequences on NCBI that are assigned to `Animalia`. 

```bash
cp /nesi/project/nesi02659/obss_2021/resources/day4/ncbi_16S.fasta ~/obss_2021/edna/references
```

Let's look at the top 10 lines of the sequencing file to determine the format.

```bash
head -10 ../references/ncbi_16S.fasta
```

```
>EF646645
ATTTAGCTAGTAACAGCAAGCAAAACGCAATTTTAGTTTGCACCCCCGAAACTAAGTGAGCTACTTTAAAACAGCCATATGGGCCAACTCGTCTCTGTTGCAAAAGAGTGAGAAGATTTTTAAGTAGAAGTGACAAGCCTATCGAACTTAGAGATAGCTGGTTATTTGGGAAAAGAATATTAGTTCTGCCTTAAGCTTTTATTAACACCCTTCAAAGCAACTTAAGAGTTATTCAAATAAGGTACAGCCTATTTGATACAGGAAACAACCTAAAACATAGGGTAACGCTACCTACAATTTTTATAATTAAGTTGGCCTAAAAGCAGCCATCTTTTAAAAAGCGTCAAAGCTTAATTATTTAATAACAACAATCACAACAAAAATGCCAAACCCACCACTACTACCAAATAACTTTATATAATTATAAAAGATATTATGCTAAAACTAGTAATAAGAAAATGACTTTCTCCTAAATATAAGTGTAATGCAGAATAAACAAATCACTGCTAATTATTGTTTTTGATTAAATAGTAGCAACCTCTCTAGAAAACCCTACTATAAATAACATTAACCTAACACAAGAATATTACAGGAAAAATTAAAAG
>EF646644
ATTTAGCTAGTCAAAACAAGCAAAGCGTAACTTAAGTTTGCTTTCCCGAAACTAAGTGAGCTACTTTGAAACAGCCTTACGGGCCAACTCGTCTCTGTTGCAAAAGAGTGAGAAGATTTTTAAGTAGAAGTGAAAAGCCTATCGAACTTAGAGATAGCTGGTTATTTGGGAAAAGAATATTAGTTCTGCCTTAAGCCAAACAACACAGTTTAAAGTAACTTAAGAGTTATTCAAATAAGGTACAGCCTATTTGATACAGGAAACAACCTAGAGAGTAGGGTAACGTCATTTAAAACGTAATTAAGTTGGCCTAAAAGCAGCCATCTTTTAAAAAGCGTCAAAGCTTAATTATTTAATAATACTAATTTCAACAAAAATACAAAACCCACATTTAATACCAAATAACTTTATATAATTATAAAAGATATTATGCTAAAACTAGTAATAAGAAAAAGACTTTCTCCTAAATATAAGTGTAACGCAGAATAAACAAATTACTGCTAATTACTGTTCATGACTCAATAGTAGTAACCTCACTAGAAAACCCTACTAATTATAACATTAATCTAACACAAGAGTATTAC
```
{: .output}


> ## Step 2: Extracting amplicon regions through *in silico* PCR analysis
>Once we have downloaded the sequencing data, we can start with a first curation step of the reference database. We will be extracting the amplicon region from each sequence through an *in silico* PCR analysis. With this analysis, we will locate the forward and reverse primer binding regions and extract the sequence in between. This will significantly reduce file sizes when larger data sets are downloaded, while also keeping only the necessary information.
>
>When we call the help documentation of the `insilico_pcr` function, we can see what parameters we need to fill in.
>
>```bash
>crabs_v1.0.1 insilico_pcr -h
>```
>
>`--fwd`: the forward primer sequence --> GACCCTATGGAGCTTTAGAC
>
>`--rev`: the reverse primer sequence --> CGCTGTTATCCCTADRGTAACT
>
>`--input`: input filename
>
>`--output`: output filename
>
>To run this command, let's create a new script using the `nano` text editor. Remember to source the `eDNA.sh` file and direct the script to the correct folder.
>
>```bash
>nano insilico.sh
>```
>
>We can exit out of `nano` by pressing `ctr+x`, followed by `y` and `enter` to save the file. Remember that to run the script, we use `bash`, followed by the filename.
>
>```bash
>bash insilico.sh
>```
>
>> ## Solution
>> Your script should look like this:
>>
>> ~~~
>> source eDNA.sh
>>
>> cd ../references/
>>
>> crabs_v1.0.1 insilico_pcr -f GACCCTATGGAGCTTTAGAC \
>> -r CGCTGTTATCCCTADRGTAACT \
>> -i ncbi_16S.fasta \
>> -o ncbi_16S_insilico.fasta
>> ~~~
>>
> {: .solution}
{: .challenge}

To check if everything executed properly, we can list the files that are in the `references` folder. We can also check how the `insilico_pcr` function has altered our file and count how many sequences we have within `ncbi_16S_insilico.fasta`.

```bash
ls -ltr ../references/
head ../references/ncbi_16S_insilico.fasta
grep -c "^>" ../references/ncbi_16S_insilico.fasta
```

> ## Step 3: Assigning a taxonomic lineage to each sequence
>Before we continue curating the reference database, we will need assign a taxonomic lineage to each sequence. The reason for this is that we have the option to curate the reference database on the taxonomic ID of each sequence. For this module, we normally need to download the taxonomy files from NCBI using the `db_download --source taxonomy` function. However, these files are already available to you for this tutorial. Therefore, we can skip the download of these files. CRABS extracts the necessary taxonomy information from these NCBI files to provide a taxonomic lineage for each sequence using the `assign_tax` function.
>
>Let's call the help documentation of the `assign_tax` function to see what parameters we need to fill in.
>
>```bash
>crabs_v1.0.1 assign_tax -h
>```
>
>`--input`: input filename 
>
>`--output`: output filename
>
>`--acc2tax`: NCBI downloaded file containing accession and taxonomic ID combinations --> nucl_gb.accession2taxid
>
>`--taxid`: NCBI downloaded file containing taxonomic ID of current and parent taxonomic level --> nodes.dmp
>
>`--name`: NCBI downloaded file containing taxonomic ID and phylogenetic name combinations --> names.dmp
>
>To run this command, let's create a new script using the `nano` text editor. Remember to source the `eDNA.sh` file and direct the script to the correct folder.
>
>```bash
>nano assign_tax.sh
>```
>
>We can exit out of `nano` by pressing `ctr+x`, followed by `y` and `enter` to save the file. Remember that to run the script, we use `bash`, followed by the filename.
>
>```bash
>bash assign_tax.sh
>```
>
>> ## Solution
>> Your script should look like this:
>>
>> ~~~
>> source eDNA.sh
>>
>> cd ../references/
>>
>> crabs_v1.0.1 assign_tax \
>>-i ncbi_16S_insilico.fasta \
>>-o ncbi_16S_insilico_tax.fasta \
>>-a nucl_gb.accession2taxid \
>>-t nodes.dmp \
>>-n names.dmp
>> ~~~
>>
> {: .solution}
{: .challenge}

We can conduct the same three checks as we did before to determine if everything executed properly. Note that we cannot use the `grep` code anymore, since the format of the file has changed. Instead, we will be counting the number of lines using the `wc -l` function.

```bash
ls -ltr ../references/
head ../references/ncbi_16S_insilico_tax.fasta
wc -l ../references/ncbi_16S_insilico_tax.fasta
```

> ## Step 4: Curating the reference database
>Next, we will run two different functions to curate the reference database, the first is to dereplicate the data using the `dereplicate` function, the second is to filter out sequences that we would not want to include in our final database using the `seq_cleanup` function. Try yourself to create a script that will run both functions with the hints and tips below.
>
> - We will call the script `curate.sh`
> - We will chose the `uniq_species` method of dereplication
> - We will use the default values within `seq_cleanup`
> - We will discard `environmental` sequences, as well as sequences that do not have a `species` name, and sequences with missing `taxonomic information`.
>
>> ## Solution
>> Your script should look like this:
>>
>> ~~~
>> source eDNA.sh
>>
>> cd ../references/
>>
>> crabs_v1.0.1 dereplicate \
>>-i ncbi_16S_insilico_tax.fasta \
>>-o ncbi_16S_insilico_tax_derep.fasta \
>>-m uniq_species
>>
>> crabs_v1.0.1 seq_cleanup \
>>-i ncbi_16S_insilico_tax_derep.fasta \
>>-o ncbi_16S_insilico_tax_derep_clean.fasta \
>>-e yes \
>>-s yes \
>>-na 0
>> ~~~
>>
> {: .solution}
{: .challenge}

We can conduct the same three checks as we did before to determine if everything executed properly.

```bash
ls -ltr ../references/
head ../references/ncbi_16S_insilico_tax_derep_clean.fasta
wc -l ../references/ncbi_16S_insilico_tax_derep_clean.fasta
```

> ## Step 5: Exporting the reference database
>Lastly, we can export the reference database to one of several formats that are implemented in the most commonly used taxonomic assignment tools using the `tax_format` function. As you have used the sintax method previously, we will export the reference database in this format. Try yourself to create a script that will export the database we created to sintax format.
>
>> ## Solution
>> Your script should look like this:
>>
>> ~~~
>> source eDNA.sh
>>
>> cd ../references/
>>
>> crabs_v1.0.1 tax_format \
>>-i ncbi_16S_insilico_tax_derep_clean.fasta \
>>-o ncbi_16S_insilico_tax_derep_clean_sintax.fasta \
>>-f sintax
>> ~~~
>>
> {: .solution}
{: .challenge}

Let's open the reference database to see what the final database looks like!

```bash
head ../references/ncbi_16S_insilico_tax_derep_clean_sintax.fasta
```

> ## Optional extra 1
>
>Now that we have created our own reference database for fish of the 16S rRNA gene. Let's try to run the taxonomic assignment again that we ran in the previous section, but now with the newly generated reference database. How do the results compare?
{: .challenge}

> ## Optional extra 2
>
>Try to download and create a reference database for shark sequences of the 16S rRNA gene. Once completed, you can run the taxonomic assignment of the previous section with this database as well. What are the results?
{: .challenge}

{% include links.md %}
