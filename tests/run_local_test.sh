#!/usr/bin/env bash
# Test the pipeline installation end to end on the bundled test data.
# Run from anywhere; requires snakemake >= 9 and conda/mamba on PATH.
#   bash tests/run_local_test.sh            # full test; tests/results/ is removed if it passes
#   bash tests/run_local_test.sh --keep     # full test, keep tests/results/ afterwards
#   bash tests/run_local_test.sh --dry-run  # only lint + DAG check
# On failure tests/results/ is always kept for debugging.
set -euo pipefail

DRY_RUN=false
KEEP=false
for arg in "$@"; do
	case "$arg" in
		--dry-run) DRY_RUN=true ;;
		--keep) KEEP=true ;;
		*) echo "Unknown option: $arg" >&2; exit 2 ;;
	esac
done

cd "$(dirname "$0")/.."
PROFILE=tests/profile
OUT=tests/results
N_SAMPLES=$(grep -vc '^sample' tests/config/metadata.tsv)

step() { printf '\n=== %s ===\n' "$*"; }
fail() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }

step "Checking tools"
command -v snakemake >/dev/null || fail "snakemake not found on PATH"
command -v conda >/dev/null || command -v mamba >/dev/null || fail "conda/mamba not found on PATH"
snakemake --version

step "Lint (warnings only)"
snakemake --profile "$PROFILE" --lint || echo "(lint reported issues; continuing)"

step "Dry run"
snakemake --profile "$PROFILE" -n --quiet rules

$DRY_RUN && { echo "Dry run OK"; exit 0; }

step "Removing old test results"
rm -rf "$OUT"

step "Creating conda environments"
snakemake --profile "$PROFILE" --conda-create-envs-only

step "Full run"
snakemake --profile "$PROFILE"

step "Checking outputs"
matrix="$OUT/counts/count_matrix.tsv"
[[ -s "$matrix" ]] || fail "$matrix missing or empty"
ncol=$(head -1 "$matrix" | awk -F'\t' '{print NF}')
(( ncol == N_SAMPLES + 1 )) || fail "$matrix has $ncol columns, expected $((N_SAMPLES + 1))"
total=$(awk -F'\t' 'NR>1{for(i=2;i<=NF;i++)s+=$i} END{print s+0}' "$matrix")
(( total > 0 )) || fail "$matrix has zero total counts"
echo "count matrix: $(($(wc -l < "$matrix") - 1)) genes, $total assigned reads"

nkraken=$(ls "$OUT"/kraken/*.kraken.out 2>/dev/null | wc -l)
(( nkraken == N_SAMPLES )) || fail "found $nkraken kraken reports, expected $N_SAMPLES"

[[ -s "$OUT/qc/multiqc.html" ]] || fail "$OUT/qc/multiqc.html missing"

step "All checks passed"

if $KEEP; then
	echo "Results kept in $OUT/"
else
	rm -rf "$OUT"
	echo "Removed $OUT/ (use --keep to keep it)"
fi
