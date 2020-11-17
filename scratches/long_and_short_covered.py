#  python long_and_short_covered.py short_reads_abundance.tsv  long_reads_abundance.tsv 1.0 50.0 isoforms_chr*.fa

import sys

from pathlib import Path

from Bio import SeqIO


def get_tpm_ind(tsv):
    with open(tsv, 'r') as fin:
        columns = fin.readline().strip().split()
    for ind, name in enumerate(columns):
        if 'tpm' in name or 'TPM' in name:
            return -(len(columns) - ind)

def get_covered_t_ids(tsv, min_tpm, max_tpm):
    ids = set()
    tpm_count_ind = get_tpm_ind(tsv)
    with open(tsv, 'r') as fin:
        next(fin)
        for line in fin:
            values = line.strip().split()
            id = values[0]
            tpm = float(values[tpm_count_ind])
            if tpm > min_tpm and tpm <= max_tpm:
                ids.add(id)
    print(tsv + ': {}'.format(len(ids)))
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
max_est = float(sys.argv[4])
in_fasta = sys.argv[5]

out_fasta = Path(in_fasta).stem + '.simultaneously_covered.{}-{}.fa'.format(min_est, max_est)

short_covered = get_covered_t_ids(short_tsv, min_est, max_est)
long_covered = get_covered_t_ids(long_tsv, min_est, max_est)

together_covered = short_covered.intersection(long_covered)
print('Simultaneously covered: {}'.format(len(together_covered)))

filter_fasta_by_ids(in_fasta, out_fasta, together_covered)