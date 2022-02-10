# FastQC Quality Check

# For each lane, pool all of the forward and reverse reads together and then run them on FastQC

rule fastqc_lane:
	input:
	output:
		html = "results/QC/lane" 
		zip = 

	conda: