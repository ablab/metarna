# python get_t_names_by_gene.py Homo_sapiens.GRCh38.82.cleared.gtf 15 listname

import sys

import random

gtf = sys.argv[1]
t_needed = int(sys.argv[2])


gene_set = set()
transcript_set = set()
chr_set = set()

lines = open(gtf).read().splitlines()
n_line = random.randint(50, len(lines) - t_needed * 30)
for line in lines[n_line:]:
    fields = line.strip().split('\t')
    chr = fields[0]
    type = fields[2]
    others = fields[8].split('; ')
    g_id = others[0].split('"')[1]
    t_id = others[2].split('"')[1]
    if type == 'transcript':
        t_needed -= 1
        gene_set.add(g_id)
        transcript_set.add(t_id)
        chr_set.add(chr)
    if t_needed == 0:
        break

print('Line: {} / {}'.format(n_line, len(lines)))
print('Genes: ' + ' '.join(gene_set))
print('Transcripts: ' + ' '.join(transcript_set))
print('chr: ' + ' '.join(chr_set))

fout = open(sys.argv[3] + '.names', 'w')
fout.write('\n'.join(transcript_set))
fout.close()