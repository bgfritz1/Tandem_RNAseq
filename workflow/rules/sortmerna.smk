
rule get_rRNA_db:
	output: 
		"resources/rRNA_ref/{ref_db}.fastq"
	log: 
		"logs/sortmerna/ref_download/{ref_db}_download.log"
	shell:
		""" curl 'https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/{ref_db}.fasta' > resources/rRNA_ref/{ref_db}.fastq 2>{log} """

rule remove_rRNA:
	input:
		fastq1="results/trimmed/{sample}-{unit}-R1.fastq.gz",
		fastq2="results/trimmed/{sample}-{unit}-R2.fastq.gz"
	output:
		fastq1="results/rRNA_dep/{sample}-{unit}_rRNAdep_fwd.fq.gz",
		fastq2="results/rRNA_dep/{sample}-{unit}_rRNAdep_rev.fq.gz",
		wd = temp(directory("results/rRNA_dep/{sample}-{unit}-temp"))
	conda:
		"../envs/sortmerna.yaml"
	log:
		"results/logs/rRNA_dep/{sample}-{unit}.log"
	params:
		workdir="results/rRNA_dep/{sample}-{unit}-temp",
		outdir="results/rRNA_dep/{sample}-{unit}_rRNAdep",
		refs=get_rNAseq_database_arg()
	shell:
		"""
		sortmerna --workdir {params.workdir} {params.refs}  \
		--reads {input.fastq1} --reads {input.fastq2} --fastx --other {params.outdir} --out2 --paired_out  1> {log}

		mv ./{params.workdir}/out/aligned.log ./{params.outdir}.log
		"""




