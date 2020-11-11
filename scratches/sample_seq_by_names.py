# sample_isoforms_by_chr.py isoforms_fa gtf chr output

from Bio import SeqIO

import sys

fasta = sys.argv[1]
names_txt = sys.argv[2]
filtered_fasta = sys.argv[3]

with open(names_txt, 'r') as fin:
    names = fin.read().splitlines()

filtered = []
for record in SeqIO.parse(fasta, "fasta"):
    if record.id[:-11] in names:
        filtered.append(record)

SeqIO.write(filtered, filtered_fasta, "fasta")