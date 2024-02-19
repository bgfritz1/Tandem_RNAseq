import sys
import pandas as pd


# logging
sys.stderr = open(snakemake.log[0], "w")


# def get_column(strandedness):
#     if pd.isnull(strandedness) or strandedness == "none":
#         return 1  # non stranded protocol
#     elif strandedness == "yes":
#         return 2  # 3rd column
#     elif strandedness == "reverse":
#         return 3  # 4th column, usually for Illumina truseq
#     else:
#         raise ValueError(
#             (
#                 "'strandedness' column should be empty or have the "
#                 "value 'none', 'yes' or 'reverse', instead has the "
#                 "value {}"
#             ).format(repr(strandedness))
#         )


counts = [
    pd.read_table(
        f, index_col=0, usecols=[0, 6], header=None, skiprows=[0,1]
    )
    for f in snakemake.input
]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")