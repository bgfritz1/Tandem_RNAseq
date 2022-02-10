rule cutadapt_pe:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}-{unit}-R1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}-R2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        "{}".format(config["params"]["cutadapt-pe"]),
    log:
        "results/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "0.31.1/bio/cutadapt/pe"