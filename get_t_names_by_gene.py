# python get_t_names_by_gene.py Homo_sapiens.GRCh38.82.cleared.gtf ENSG00000230021

import sys

gtf = sys.argv[1]
gene = sys.argv[2]

fout = open(gene + '.names', 'w')

with open(gtf, 'r') as fin:
    for line in fin:
        fields = line.strip().split('\t')
        type = fields[2]
        others = fields[8].split('; ')
        g_id = others[0].split('"')[1]
        t_id = others[2].split('"')[1]
        if type == 'transcript' and g_id == gene:
            fout.write(t_id + '\n')

fout.close()