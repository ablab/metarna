# sample_isoforms_by_chr.py isoforms_fa gtf chr output

from Bio import SeqIO
from BCBio import GFF

import sys

fasta = sys.argv[1]
gtf = sys.argv[2]
chr = sys.argv[3]
filtered_fasta = sys.argv[4]

# sample isoform ids which is chr owned from gtf file
limit_info = dict(gff_id=[chr], gff_type = ["transcript"])
in_handle = open(gtf)
ids = []
for rec in GFF.parse(in_handle, limit_info=limit_info):
    for transcript in rec.features:
        ids.append(transcript.qualifiers['transcript_id'][0])
in_handle.close()

# choose this isoform from fasta file
filtered = []
for record in SeqIO.parse(fasta, "fasta"):
    if record.id[:-11] in ids:
        filtered.append(record)

SeqIO.write(filtered, filtered_fasta, "fasta")