import os
import subprocess

path = '/Nancy/data/input/RNA/ENA_gut/alex/'
bar = '/Nancy/data/input/RNA/ENA_gut/alex/bar.txt'
outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

all_proteins = []
with open(bar, 'r') as fin:
    lines = fin.readlines()
    for mode in ['rna', 'meta', 'rnaviral']:
        command = 'cat'
        all_proteins.append(os.path.join(outdir, 'all.{mode}.proteins.faa'.format(mode=mode)))
        for i in range(0, len(lines) - 1, 2):
            name = lines[i].strip().split('/')[-2]
            proteins = os.path.join(outdir, 'prodigal_out', '{name}.{mode}.proteins.faa'.format(name=name, mode=mode))
            command += ' ' + proteins
        command += ' > {all_proteins}'.format(all_proteins=all_proteins[-1])
        #print(command)
        subprocess.call(command, shell=True)
    #print('\n')

    for proteins in all_proteins:
        name = proteins.split('.')[1]
        #print(name)
        command = 'mmseqs easy-linclust {all_proteins} {clusterRes} {tmp} ' \
                  '--min-seq-id 0.5 -c 0.5 --cov-mode 1 --cluster-mode 2 --kmer-per-seq 80'.\
            format(all_proteins=proteins, clusterRes=name, tmp=os.path.join(outdir, 'tmp'))
        #print(command)
        subprocess.call(command, shell=True)
    #print('\n')

