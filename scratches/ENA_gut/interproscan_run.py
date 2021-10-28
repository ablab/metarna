import os
import subprocess

#alex_dir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'
ipr_dir = '/Nancy/ebushmanova/ENA_gut/alex_runs/interproscan_out_1/'

for mode in ['rna', 'meta', 'rnaviral']:
    cat_cmd = 'cat'
    mode_dir = os.path.join(ipr_dir, '{}_rep_seq.clear'.format(mode))
    for num in ['0' + str(n) for n in range(0, 10)] + list(range(10, 50)):
        #print(num)
        proteins = os.path.join(mode_dir, '{mode}_rep_seq.clear.{num}.fasta'.format(mode=mode, num=num))
        #command = 'interproscan.sh -i {proteins} -d {dir} -cpu 16 -dp'.format(proteins=proteins, dir=mode_dir)
        #subprocess.call(command, shell=True)
        #print(command)
        cat_cmd += ' {}.tsv'.format(proteins)
    cat_cmd += ' > ' + os.path.join(ipr_dir, '{}_rep_seq.clear.tsv'.format(mode))
    subprocess.call(cat_cmd, shell=True)
    #print(cat_cmd)
