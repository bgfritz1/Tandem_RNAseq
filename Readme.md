# README

## To run the workflow, do the following: 

1. Make sure that the necessary prerequisite files are downloaded. At least this should be the human genome fastq and the annotation file
  - FASTQ: http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  - Annotation: http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.chr.gtf.gz

2. Make sure the kraken database is downloaded and accessible. Prebuilt databases can be found here: https://benlangmead.github.io/aws-indexes/k2

3. Edit the `metadata.tsv` file with your sample information. The metadata file is tab-delimited matrix and should at least have the following column values:
  - `sample`: Where every sample should have a unique sample name
  - `strategy`: either `PAIRED` or `SINGLE` referring the whether it was paired-end sequencing or single-end. 
  - `strandedness`: refers to the a stranded or directional library prep. `0 = unstranded, 1 = forward-stranded, 2=reverse-stranded`. 
  - Additonal columns are also accepted.
 
 4. Edit the `units.tsv` file. This contains information about the raw sequence data files. Paths can be either relative or absolute. The `unit`column refers
 to samples that were sequenced across multiple lanes on the same flow cell. If the reads are single end, just leave the `fq2` column blank.
 
 5. Edit the `config/congif.yaml`file with the relevant details. Paths can be relative to the repository directory. Note: Some options in the config file
 have not yet been implemented.
 
 6. Once everythng is set up, make sure the following programs are available:
  - `snakemake v6.15.3` I have had problems running the pipeling on snakemake v7
  - `mamba 0.24.0`
  - `python>3.7`  
  On computerome, this can be done with: 
  ```
  module load miniconda3/4.12.0
  module load mamba-org/mamba/0.24.0
  module load snakemake/6.15.3
  ```
  
  7. In `config/Computerome2/config.yaml` change the `cluster` parameter to the proper Computerome account string. 
  
  8. Now, you should be read to go! Navigate to the repository directory (the one containing the config/, workflow/, etc. files) and run the following command: 
  `snakemake --profile config/Computerome2`. This will start the analysis. It is advisable to run the analysis on a subset of samples first 
  to check that it works as expected. 
