import os

import subprocess

from Bio import SeqIO

outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs'

def get_annotated_proteins(align_path, unique):
    proteins = set()
    with open(align_path, 'r') as fin:
        for line in fin:
            str = line.strip().split()[0]
            if unique:
                ind = str.rfind('_')
                proteins.add(str[:ind])
            else:
                proteins.add(str)
    return proteins

def get_annotation_stats(outdir, suffix):
    print(suffix + '\n')
    for mode in ['rna', 'meta', 'rnaviral']:
        path = os.path.join(outdir, mode + suffix)
        annotated = len(get_annotated_proteins(path, False))
        annotated_ignoring = len(get_annotated_proteins(path, True))
        print('{}\n#annotated: {}\n#annotated ignoring _: {}\n\n'.format(mode, annotated, annotated_ignoring))

get_annotation_stats(os.path.join(outdir, 'interproscan_out_1'), '_rep_seq.clear.tsv')
get_annotation_stats(outdir, '.matches.m8')
