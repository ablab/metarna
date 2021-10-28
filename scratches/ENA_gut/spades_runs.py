import os
import subprocess

path = '/Nancy/data/input/RNA/ENA_gut/alex/'
bar = '/Nancy/data/input/RNA/ENA_gut/alex/bar.txt'
spades = '/home/ebushmanova/algorithmic-biology/assembler/spades.py'
outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

with open(bar, 'r') as fin:
    lines = fin.readlines()
    for i in range(0, len(lines) - 1, 2):
        #print(i)
        name = lines[i].strip().split('/')[-2]
        left = os.path.join(path, lines[i].strip().split('/')[-1])
        right = os.path.join(path, lines[i + 1].strip().split('/')[-1])

        for mode in ['rna', 'meta', 'rnaviral']:
            command = '{spades} --{mode} -1 {left} -2 {right} -o {outdir} -t 64'.\
                format(spades=spades, mode=mode, left=left, right=right,
                       outdir=os.path.join(outdir, mode + 'spades_' + name))
            subprocess.call(command, shell=True)
            #print(command)
        #print('/n')

