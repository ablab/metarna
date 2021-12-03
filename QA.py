#!/usr/bin/env python

# Usage:
# ./QA.py path/to/transcripts/ path/to/output/directory

import os
import subprocess
import sys

from Bio import SeqIO


threads = 16

mgy_db = '/Nancy/data/input/RNA/ENA_gut/db/mgy.dmnd'


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


def run_prodigal(assembly_path, outdir):
    initial_dir = os.path.abspath(os.getcwd())
    os.chdir(outdir)
    name = os.path.basename(assembly_path).split('.')[0]
    proteins = '{}.proteins.faa'.format(name)
    command = 'prodigal -i {assembly} -o {genes} -a {proteins} -p meta'.\
        format(assembly=assembly_path, genes='{}.genes.faa'.format(name), proteins=proteins)
    print(command)
    subprocess.call(command, shell=True)
    os.chdir(initial_dir)
    return os.path.join(outdir, proteins)


def run_mmseqs(proteins_path, outdir, min_seq_id=0.9):
    name = os.path.basename(proteins_path).split('.')[0]
    rep_seq_path = os.path.join(outdir, '{}_rep_seq.fasta'.format(name))
    command = 'mmseqs easy-linclust {proteins} {clusterRes} {tmp} ' \
              '--min-seq-id {min_seq_id} --cov-mode 1 --cluster-mode 2 --kmer-per-seq 80'. \
        format(proteins=proteins_path, clusterRes=os.path.join(outdir, name), min_seq_id=min_seq_id, tmp=os.path.join(outdir, 'tmp'))
    print(command)
    subprocess.call(command, shell=True)
    return rep_seq_path


def run_interproscan(rep_seq_path, outdir):
    global threads

    name = os.path.basename(rep_seq_path).split('_')[0]

    ipr_dir = os.path.join(outdir, '{}_rep_seq.clear'.format(name))
    if not os.path.exists(ipr_dir):
        os.makedirs(ipr_dir)

    # Remove the * at the end of the sequences
    clear_rep_seq = os.path.join(ipr_dir, '{}_rep_seq.clear.fasta'.format(name))
    command = 'sed \'s/*//\' ' + rep_seq_path + ' > ' + clear_rep_seq
    print(command)
    subprocess.call(command, shell=True)

    # Split fasta into smaller chunks to speed up
    command = 'pyfasta split -n 50 {}'.format(clear_rep_seq)
    print(command)
    subprocess.call(command, shell=True)

    cat_cmd = 'cat'
    for num in ['0' + str(n) for n in range(0, 10)] + list(range(10, 50)):
        # print(num)
        filename = '{name}_rep_seq.clear.{num}'.format(name=name, num=num)
        proteins_path = '{}.fasta'.format(filename)
        # Run IPR for each chunk separately
        command = 'interproscan.sh -i {proteins} -b {base} -cpu {threads} -dp'.format(proteins=os.path.join(ipr_dir, proteins_path), base=os.path.join(ipr_dir, filename), threads=threads)
        print(command)
        subprocess.call(command, shell=True)
        cat_cmd += ' {}.tsv'.format(os.path.join(outdir, proteins_path))
    ipr_path = os.path.join(ipr_dir, '{}_rep_seq.clear.tsv'.format(name))
    # Concatenate all IPR tsv-s
    cat_cmd += ' > ' + ipr_path
    print(cat_cmd)
    subprocess.call(cat_cmd, shell=True)
    return ipr_path


def run_diamond(rep_seq_path, mgy_db, outdir):
    global threads

    name = os.path.basename(rep_seq_path).split('_')[0]
    diamond_path = os.path.join(outdir, '{}.matches.m8'.format(name))
    command = 'diamond blastp -d {mgy} -q {rep_seq} -o {diamond} --query-cover 50 --id 95 --subject-cover 90 ' \
              '--threads {threads} -dra -appl Hamap,Pfam'.\
        format(mgy=mgy_db, rep_seq=rep_seq_path, diamond=diamond_path, name=name, threads=threads)
    print(command)
    subprocess.call(command, shell=True)
    return diamond_path


def get_counts(rep_seq_path, diamond_path, ipr_path):
    annotated_mgy = get_annotated_proteins(diamond_path)
    annotated_ipr = get_annotated_proteins(ipr_path)
    all_clusters = get_fasta_ids(rep_seq_path)
    # print(list(annotated_ipr)[:5], list(annotated_mgy)[:5], list(all_clusters)[:5])
    none_proteins = all_clusters - annotated_mgy - annotated_ipr
    both_proteins = annotated_mgy.intersection(annotated_ipr)
    mgy_only = annotated_mgy - both_proteins
    ipr_only = annotated_ipr - both_proteins
    sum = len(none_proteins) + len(mgy_only) + len(ipr_only) + len(both_proteins)
    print('MGnify: {}\nIPR: {}'.format(len(annotated_mgy), len(annotated_ipr)))
    print('None: {}\nMGnify only: {}\nIPR only: {}\nboth: {}\nsummary: {}'
          .format(len(none_proteins), len(mgy_only), len(ipr_only), len(both_proteins), sum))
    print('All clusters: {}\n'.format(len(all_clusters)))


def main():
    global mgy_db

    assembly_path = sys.argv[1]
    outdir = sys.argv[2]

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    proteins_path = run_prodigal(assembly_path, outdir)

    rep_seq_path = run_mmseqs(proteins_path, outdir)

    ipr_path = run_interproscan(rep_seq_path, outdir)

    diamond_path = run_diamond(rep_seq_path, mgy_db, outdir)

    get_counts(rep_seq_path, diamond_path, ipr_path)


if __name__ == '__main__':
    main()
