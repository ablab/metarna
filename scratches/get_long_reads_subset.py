#python get_long_reads_subset.py rand25.names Homo_sapiens.GRCh38.82.cleared.isoforms.fa /Nancy/ebushmanova/mt_feb/long_reads_data/hybridRNA/MCF7/isoseq_all.fastq.gz

import sys

import subprocess

from pathlib import Path

commands = []

names_txt = sys.argv[1]
isoforms_fasta = sys.argv[2]
long_reads = sys.argv[3]

names_filename = Path(names_txt).stem.split('.')[0]
long_reads_filename = Path(long_reads).stem.split('.')[0]


commands.append('minimap2 -x map-pb -t 16 -a -L {} {} > {}.Aligned.sam'.format(isoforms_fasta, long_reads, long_reads_filename))
commands.append('samtools view -b {name}.Aligned.sam > {name}.Aligned.bam'.format(name=long_reads_filename))
commands.append('samtools sort {name}.Aligned.bam > {name}.Aligned.sortedByCoord.bam'.format(name=long_reads_filename))
commands.append('samtools index {}.Aligned.sortedByCoord.bam'.format(long_reads_filename))

with open(names_txt, 'r') as fin:
    names = fin.read().splitlines()
print('seqnames: {}\n'.format(names))

for seqname in names:
    commands.append('samtools view -h -b {name}.Aligned.sortedByCoord.bam {seqname}_transcript > {seqname}.bam'.
                    format(name=long_reads_filename, seqname=seqname))
    commands.append('samtools sort -n {seqname}.bam -o {seqname}.sortedByName.bam'.format(seqname=seqname))
    commands.append('samtools bam2fq {seqname}.sortedByName.bam > {name}.{seqname}.fastq'.
                    format(name=long_reads_filename, seqname=seqname))

reads = ['{}.{}.fastq'.format(long_reads_filename, seqname) for seqname in names]
commands.append('cat {} > {}.{}.fastq'.format(' '.join(reads), long_reads_filename, names_filename))

for cmd in commands:
    print(cmd + ':\n')
    subprocess.call(cmd, shell=True)
    print('\n')