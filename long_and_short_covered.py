#  python long_and_short_covered.py short_reads_abundance.tsv  long_reads_abundance.tsv 50.0 isoforms_chr*.fa

import sys

from pathlib import Path

from Bio import SeqIO


def get_covered_t_ids(tsv, est_count_idx, min_est):
    ids = set()
    with open(tsv, 'r') as fin:
        next(fin)
        for line in fin:
            values = line.strip().split()
            id = values[0]
            est = float(values[est_count_idx])
            if est >= min_est:
                ids.add(id)
    return ids

def filter_fasta_by_ids(in_fasta, out_fasta, ids):
    records = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in ids:
            records.append(record)
    SeqIO.write(records, out_fasta, "fasta")
    return out_fasta


short_tsv = sys.argv[1]
long_tsv = sys.argv[2]
min_est = float(sys.argv[3])
in_fasta = sys.argv[4]

out_fasta = Path(in_fasta).stem + '.simultaneously_covered.fa'

short_covered = get_covered_t_ids(short_tsv, 3, min_est)
long_covered = get_covered_t_ids(long_tsv, 1, min_est)
together_covered = short_covered.intersection(long_covered)

filter_fasta_by_ids(in_fasta, out_fasta, together_covered)