import pandas as pd
import pdb 

#Read in the sample and unit information

samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False
)
samples.index.names = ["sample_id"]

def drop_unique_cols(df):
    singular_cols = df.nunique().loc[(df.nunique().values <= 1)].index
    return df.drop(singular_cols, axis=1)

samples = drop_unique_cols(samples)

units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(
    ["sample", "unit"], drop=False
)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
) 

##### wildcard constraints #####

wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"])

#Helper functions 

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    return [f"{u.fq1}", f"{u.fq2}"]

def get_trimmed(wildcards):
        return expand(
            "results/rRNA_dep/{sample}-{unit}_rRNAdep_{group}.fq.gz",
            group=['fwd', 'rev'],
            **wildcards,
        )

def get_rNAseq_database_arg():
    databases = config["refs"]["rRNA_ref_dbs"]
    arg = ""
    for database in databases:
        arg = f"{arg} --ref resources/rRNA_ref/{database}.fastq"
    return(arg)

def all_input(wildcards):
    """
    Function defining all requested inputs for the rule all (below).
    """
    wanted_input = []
    #Probably Not Necessary later on since something will use the cutadapt output
    if config["modules"]["cutadapt"] == "enable":
        wanted_input.extend(
            expand(
                "results/generated_data/trimmed/{unit.sample}-{unit.unit}-R{read}.fastq.gz",
                unit = units[["sample", "unit"]].itertuples(),
                read = [1,2]
            )
        )
    # Kallisto output
    wanted_input.extend(
        expand(
            "results/kallisto/{unit.sample}-{unit.unit}", 
            unit = units[["sample", "unit"]].itertuples()
            )
    )
      
    #TEST rRNA ref db
    wanted_input.extend(
        expand("resources/rRNA_ref/{ref_db}.fastq",
            ref_db = config["refs"]["rRNA_ref_dbs"]
        )
    )
