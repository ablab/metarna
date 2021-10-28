import os

import subprocess

from Bio import SeqIO

outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs'

def get_fasta_ids(fasta, unique):
    ids = set()
    for seq_record in SeqIO.parse(fasta, "fasta"):
        str = seq_record.id
        if unique:
            ind = str.rfind('_')
            ids.add(seq_record.id[:ind])
        else:
            ids.add(seq_record.id)
    return ids

for mode in ['rna', 'meta', 'rnaviral']:
    rep_seq_path = os.path.join(outdir, '{}_rep_seq.fasta'.format(mode))
    all_seqs_path = os.path.join(outdir, '{}_all_seqs.fasta'.format(mode))
    rep_seq = len(get_fasta_ids(rep_seq_path, False))
    rep_seq_ignoring = len(get_fasta_ids(rep_seq_path, True))
    all_seqs = len(get_fasta_ids(all_seqs_path, False))
    all_seqs_ignoring = len(get_fasta_ids(all_seqs_path, True))
    print('{}\n#rep_seq: {}\n#rep_seq ignoring _: {}\n'
          '#all_seqs: {}\n#all_seqs ignoring _: {}\n\n'.
          format(mode, rep_seq, rep_seq_ignoring, all_seqs, all_seqs_ignoring))
