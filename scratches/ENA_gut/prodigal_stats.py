import os

import subprocess

from Bio import SeqIO


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

bar = '/Nancy/data/input/RNA/ENA_gut/alex/bar.txt'
outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

mode_to_name = {'rna': 'transcripts', 'meta': 'scaffolds', 'rnaviral': 'scaffolds'}

with open(bar, 'r') as fin:
    lines = fin.readlines()
    for mode in ['rna', 'meta', 'rnaviral']:
        command = 'cat'
        for i in range(0, len(lines) - 1, 2):
            name = lines[i].strip().split('/')[-2]
            command += ' ' + os.path.join(outdir, 'prodigal_out', name + '.' + mode + '.proteins.faa')
        all_seq_path = os.path.join(outdir, 'all.{}.proteins.faa'.format(mode))
        command += ' > ' + all_seq_path
        subprocess.call(command, shell=True)
        # print(command + '\n')

        num = len(get_fasta_ids(all_seq_path, False))
        num_ignoring = len(get_fasta_ids(all_seq_path, True))
        print('{}\n#proteins: {}\n#proteins ignoring _: {}\n\n'.format(mode, num, num_ignoring))
