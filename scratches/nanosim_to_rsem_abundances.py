# python nanosim_to_rsem_abundances.py nanosim_abundance.tsv rsem_abundance.tsv

import sys

# nanosim abundance file format
# target_id est_counts tpm

# rsem abundance file format
# transcript_id gene_id length effective_length expected_count TPM FPKM IsoPct
# The simulator only reads the TPM column.

fin = open(sys.argv[1], 'r')
fout = open(sys.argv[2], 'w')

next(fin)
fout.write('transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n')
for line in fin:
    values = line.strip().split()
    transcript_id = values[0]
    tpm = values[2]
    fout.write('{}\t*\t*\t*\t*\t{}\t*\t*\n'.format(transcript_id, tpm))

fin.close()
fout.close()
