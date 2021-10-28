import os
import subprocess

alex_dir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

for mode in ['rna', 'meta', 'rnaviral']:
    command = 'sed \'s/*//\' ' + \
              os.path.join(alex_dir, '{}_rep_seq.fasta'.format(mode)) + ' > ' + \
              os.path.join(alex_dir, '{}_rep_seq.clear.fasta'.format(mode))
    subprocess.call(command, shell=True)
    #print(command)

