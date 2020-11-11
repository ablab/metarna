import sys

import subprocess


commands = []

# isoforms_fasta = sys.argv[2]
# left_fq = sys.argv[3]
# right_fq = sys.argv[4]

# commands.append('bwa index {}'.format(isoforms_fasta))
# commands.append('bwa mem -t 16 {} {} {} > aln-pe.bwamem.sam'.format(isoforms_fasta, left_fq, right_fq))
# commands.append('samtools view -b aln-pe.bwamem.sam > aln-pe.bwamem.bam')
# commands.append('samtools sort aln-pe.bwamem.bam > aln-pe.bwamem.sortedByCoord.bam')
# commands.append('samtools index aln-pe.bwamem.sortedByCoord.bam')

names_txt = sys.argv[1]

filename = names_txt.split('.')[0]

with open(names_txt, 'r') as fin:
    names = fin.read().splitlines()
print('seqnames: {}\n'.format(names))

for seqname in names:
    commands.append('samtools view -h -b aln-pe.bwamem.sortedByCoord.bam {name}_transcript > {name}.bam'.format(name=seqname))
    commands.append('samtools sort -n {name}.bam -o {name}.sortedByName.bam'.format(name=seqname))
    commands.append('bedtools bamtofastq -i {name}.sortedByName.bam -fq left.{name}.fastq -fq2 right.{name}.fastq'.format(name=seqname))
    commands.append('bedtools bamtofastq -i {name}.sortedByName.bam -fq left.{name}.fastq -fq2 right.{name}.fastq'.format(name=seqname))

lefts = ['left.{}.fastq'.format(name) for name in names]
commands.append('cat {} > left.{}.fastq'.format(' '.join(lefts), filename))

rights = ['right.{}.fastq'.format(name) for name in names]
commands.append('cat {} > right.{}.fastq'.format(' '.join(rights), filename))

for cmd in commands:
    print(cmd + ':\n')
    subprocess.call(cmd, shell=True)
    print('\n')