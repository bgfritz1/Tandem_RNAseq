#Workflow for Tandem Project

1. Cutadapt 
  cutadapt-pe: "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 33 -e 0.005 --overlap 7"

2. FastQC

3. SortMeRNA

4. Kallisto 

5. Tximport to Count table

6. Kraken 


# To use

Make sure Snakemake is loaded 

Make sure CDNA fasta is in resources/....
Right now it is called Homo_sapiens.GRCh38.cdna.all.fa.gz
I need to add a cnofig a