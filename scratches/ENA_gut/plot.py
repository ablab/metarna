import os
import subprocess

from Bio import SeqIO

import numpy as np
import matplotlib.pyplot as plt

outdir = '/Nancy/ebushmanova/ENA_gut/alex_runs/'

def get_fasta_ids(fasta):
    ids = set()
    for seq_record in SeqIO.parse(fasta, "fasta"):
        ids.add(seq_record.id)
    return ids

def get_annotated_proteins(align_path):
    proteins = set()
    with open(align_path, 'r') as fin:
        for line in fin:
            proteins.add(line.strip().split()[0])
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
        all_proteins_path = os.path.join(outdir, 'all.{mode}.proteins.faa'.format(mode=mode))
        diamond_path = os.path.join(outdir, 'all.{mode}.mgy'.format(mode=mode))
        ipr_path = os.path.join(outdir, 'all.{mode}.ipr'.format(mode=mode))

        annotated_mgy = get_annotated_proteins(diamond_path)
        annotated_ipr = get_annotated_proteins(ipr_path)
        all_proteins = get_fasta_ids(all_proteins_path)
        print(annotated_ipr[:5], annotated_mgy[:5], all_proteins[:5])

        none = all_proteins - annotated_mgy - annotated_ipr
        both = annotated_mgy.intersection(annotated_ipr)

        none.append(len(none))
        mgnify.append(len(annotated_mgy))
        ipr.append(len(annotated_ipr))
        both.append(len(both))
    plot(none, mgnify, ipr, both)


if __name__ == '__main__':
    main()
