#!/usr/bin/env bash
# Rebuild tests/data/fastq/ so it only contains read pairs that map to chr22.
# Needs STAR BAMs from a previous test run (tests/results/star/pe/<sample>/pe_aligned.bam,
# e.g. from: bash tests/run_local_test.sh --keep)
# and samtools on PATH. Original FASTQs are moved to tests/data/fastq_full/.
#   bash tests/make_chr22_test_reads.sh
set -euo pipefail

cd "$(dirname "$0")/.."
BAM_DIR=tests/results/star/pe
SRC=tests/data/fastq_full
DST=tests/data/fastq

if [[ ! -d "$SRC" ]]; then
	mv "$DST" "$SRC"
	mkdir -p "$DST"
fi

tail -n +2 tests/config/units.tsv | while IFS=$'\t' read -r sample unit fq1 fq2 || [[ -n "$sample" ]]; do
	bam="$BAM_DIR/$sample/pe_aligned.bam"
	ids=$(mktemp)
	# read names of mapped primary alignments (either mate)
	samtools view -F 0x904 "$bam" | cut -f1 | sort -u > "$ids"
	for fq in "$fq1" "$fq2"; do
		in="$SRC/$(basename "$fq")"
		out="$DST/$(basename "$fq")"
		zcat "$in" | awk 'NR==FNR{keep[$1]; next}
			FNR%4==1{split(substr($1,2), a, " "); p=(a[1] in keep)}
			p' "$ids" - | gzip -n > "$out"
		echo "$out: $(( $(zcat "$out" | wc -l) / 4 )) reads"
	done
	rm "$ids"
done
