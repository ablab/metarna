import os
import subprocess

bar = '/Nancy/data/input/RNA/ENA_gut/alex/bar.txt'
outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

mode_to_name = {'rna': 'transcripts', 'meta': 'scaffolds', 'rnaviral': 'scaffolds'}

with open(bar, 'r') as fin:
    lines = fin.readlines()
    for mode in ['rna', 'meta', 'rnaviral']:
        command = 'cat'
        for i in range(0, len(lines) - 1, 2):
            name = lines[i].strip().split('/')[-2]
            command += ' ' + os.path.join(outdir, mode + 'spades_' + name, mode_to_name[mode] + '.fasta')
        all_seq_path = os.path.join(outdir, mode + 'spades.all.fasta')
        command += ' > ' + all_seq_path
        subprocess.call(command, shell=True)
        # print(command + '\n')

        command = 'quast.py ' + all_seq_path + ' -o ' + os.path.join(outdir, mode + '_quast_out')
        subprocess.call(command, shell=True)
        #print(command + '\n')
