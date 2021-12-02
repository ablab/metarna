#!/bin/bash

diamond makedb --in /Nancy/data/input/RNA/ENA_gut/db/mgy_clusters.fa -d mgy --threads 16
for mode in 'rna' 'meta' 'rnaviral'
    do
        diamond blastp -d mgy -q ${mode}_rep_seq.fasta -o ${mode}.matches.m8 --query-cover 50 --id 50 --threads 16
    done
