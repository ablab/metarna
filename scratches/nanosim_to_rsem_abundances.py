# python nanosim_to_rsem_abundances.py nanosim_abundance.tsv *.isoforms.results rsem_abundance.tsv

import sys

# rsem abundance file format
# transcript_id gene_id length effective_length expected_count TPM FPKM IsoPct
# Start from a learned *.isoforms.results file and only modify the TPM column.
learned_abundance = sys.argv[2]

def get_tpm_dict(nanosim_abundance):
    # nanosim abundance file format
    # target_id est_counts tpm
    tpm_dict = {}
    with open(nanosim_abundance, 'r') as fin:
        next(fin)
        for line in fin:
            values = line.strip().split()
            transcript_id = values[0]
            tpm = values[2]
            tpm_dict[transcript_id] = tpm
    return tpm_dict

tpm_dict = get_tpm_dict(sys.argv[1])

# tst_num = 0
fout = open(sys.argv[3], 'w')
with open(learned_abundance, 'r') as fin:
    line = fin.readline()
    fout.write(line)
    for line in fin:
        values = line.strip().split()
        transcript_id = values[0]
        if transcript_id in tpm_dict:
            values[5] = tpm_dict[transcript_id]
        else:
            values[5] = '0.0000'
            # tst_num += 1
            # print(transcript_id)
        fout.write('\t'.join(values) + '\n')
# print(tst_num)
fout.close()
