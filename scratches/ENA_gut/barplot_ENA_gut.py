import os

from Bio import SeqIO

import numpy as np
import matplotlib.pyplot as plt

# outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'
outdir = '/home/letovesnoi/work/nancy_mp/ENA_gut/alex_runs/'

def get_fasta_ids(fasta):
    ids = set()
    for seq_record in SeqIO.parse(fasta, "fasta"):
        str = seq_record.id
        ind = str.rfind('_')
        ids.add(seq_record.id[:ind])
    return ids

def get_annotated_proteins(align_path):
    proteins = set()
    with open(align_path, 'r') as fin:
        for line in fin:
            str = line.strip().split()[0]
            ind = str.rfind('_')
            proteins.add(str[:ind])
    return proteins

def plot(none, mgnify, ipr, both):
    ind = np.arange(len(none))  # the x locations for the groups
    width = 0.35  # the width of the bars: can also be len(x) sequence

    p1 = plt.bar(ind, both, width)
    p2 = plt.bar(ind, ipr, width, bottom=both)
    p3 = plt.bar(ind, mgnify, width, bottom=np.array(both) + np.array(ipr))
    p4 = plt.bar(ind, none, width, bottom=np.array(both) + np.array(ipr) + np.array(mgnify))

    plt.ylabel('Number of protein clusters')
    plt.xticks(ind, ['rna', 'meta', 'rnaviral'])
    # plt.yticks(np.arange(0, 81, 10))
    plt.legend((p4[0], p3[0], p2[0], p1[0]), ('None', 'MGnify only', 'IPR only', 'both'))

    plt.show()

# plot([1, 3, 2], [1, 2, 2], [0.5, 1.5, 1], [1.5, 2.5, 3])

def main():
    none = []
    mgnify = []
    ipr = []
    both = []
    for mode in ['rna', 'meta', 'rnaviral']:
        all_clusters_path = os.path.join(outdir, '{}_rep_seq.clear.fasta'.format(mode))
        diamond_path = os.path.join(outdir, '{}.matches.m8'.format(mode))
        ipr_path = os.path.join(outdir, 'interproscan_out_1', '{}_rep_seq.clear.tsv'.format(mode))

        annotated_mgy = get_annotated_proteins(diamond_path)
        annotated_ipr = get_annotated_proteins(ipr_path)
        # annotated_ipr = set()
        all_clusters = get_fasta_ids(all_clusters_path)
        print(list(annotated_ipr)[:5], list(annotated_mgy)[:5], list(all_clusters)[:5])

        none_proteins = all_clusters - annotated_mgy - annotated_ipr
        both_proteins = annotated_mgy.intersection(annotated_ipr)
        mgy_only = annotated_mgy - both_proteins
        ipr_only = annotated_ipr - both_proteins

        none.append(len(none_proteins))
        mgnify.append(len(mgy_only))
        ipr.append(len(ipr_only))
        both.append(len(both_proteins))

        sum = len(none_proteins) + len(mgy_only) + len(ipr_only) + len(both_proteins)
        print('None: {}\nMGnify only: {}\nIPR only: {}\nboth: {}\nsummary: {}'
              .format(len(none_proteins), len(mgy_only), len(ipr_only), len(both_proteins), sum))
        print('All clusters: {}\n'.format(len(all_clusters)))

    plot(none, mgnify, ipr, both)


if __name__ == '__main__':
    main()

