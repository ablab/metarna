# python number_of_pathes.py alignment.tsv len_hist.png

import sys

import matplotlib.pyplot as plt

cnt = 0
lens = []
with open(sys.argv[1], 'r') as fin:
    for line in fin:
        path = line.split()[6].replace(';', ',').split(',')
        if len(path) > 2:
            cnt += 1
            lens.append(len(path))

print('Number of pathes longer 2: {}'.format(cnt))

bins = range(min(lens), max(lens), 1)
plt.hist(lens, color='blue', edgecolor='black', bins=bins)
plt.xlabel('len')
plt.ylabel('# pathes')
plt.gca().set_xticks(bins)
plt.savefig(sys.argv[2])