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