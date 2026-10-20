[![DOI](https://zenodo.org/badge/759815461.svg)](https://zenodo.org/badge/latestdoi/759815461) [![Pipeline test](https://github.com/bgfritz1/Tandem_RNAseq/actions/workflows/test.yml/badge.svg)](https://github.com/bgfritz1/Tandem_RNAseq/actions/workflows/test.yml)

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
   - `outdir` sets where all pipeline outputs are written (default `results`). It can be relative to the repository directory or absolute,
   and can also be overridden at runtime, e.g. `snakemake --profile config/Computerome2 --config outdir=/home/projects/ku_00132/scratch/run1`.
 
 6. Once everythng is set up, make sure the following programs are available:
  - `snakemake >= 9` (the Snakefile enforces `min_version("9.0.0")`)
  - `conda` or `mamba` (rule environments are built with `software-deployment-method: conda`)
  - `snakemake-executor-plugin-cluster-generic` (only needed for the Computerome2 profile)

  7. In `config/Computerome2/config.yaml` change the `cluster` parameter to the proper Computerome account string. 
  
  8. Now, you should be read to go! Navigate to the repository directory (the one containing the config/, workflow/, etc. files) and run the following command: 
  `snakemake --profile config/Computerome2`. This will start the analysis. It is advisable to run the analysis on a subset of samples first 
  to check that it works as expected. 


## Testing the installation

A small test data set is bundled in `resources/` (4 paired-end samples, the human chr22 reference + annotation, and a tiny Kraken2 database).
`config/test/config.yaml` points the pipeline at it and `config/local` is a profile for running locally without a cluster.

```
mamba create -n tandem -c conda-forge -c bioconda 'snakemake>=9'
conda activate tandem
bash tests/run_local_test.sh            # lint, dry run, build envs, full run, check outputs
bash tests/run_local_test.sh --dry-run  # quick check of the DAG only
```

Outputs are written to `results_test/`. The first run takes longer because the conda environments are built.
To test on Computerome instead: `snakemake --profile config/Computerome2 --configfile config/test/config.yaml`.

The test FASTQs in `resources/test_data/` only contain read pairs that align to chr22 (they were subset with `tests/make_chr22_test_reads.sh`), so the full test runs in about 1.5 minutes once the conda environments exist.
