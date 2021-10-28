import os
import subprocess

path = '/Nancy/data/input/RNA/ENA_gut/alex/'
bar = '/Nancy/data/input/RNA/ENA_gut/alex/bar.txt'
outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

mode_to_name = {'rna': 'transcripts', 'meta': 'scaffolds', 'rnaviral': 'scaffolds'}
with open(bar, 'r') as fin:
    lines = fin.readlines()
    for i in range(0, len(lines) - 1, 2):
        name = lines[i].strip().split('/')[-2]
        for mode in mode_to_name.keys():
            transcripts = os.path.join(outdir, mode + 'spades_' + name, mode_to_name[mode] + '.fasta')
            command = 'prodigal -i {transcripts} -o {name}.{mode}.genes -a {name}.{mode}.proteins.faa -p meta'.\
                format(transcripts=transcripts, name=name, mode=mode)
            subprocess.call(command, shell=True)
            #print(command)
        #print('\n')

