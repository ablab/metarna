import sys

cnt = 0
with open(sys.argv[1], 'r') as fin:
    for line in fin:
        path = line.split()[6].replace(';', ',').split(',')
        if len(path) > 2:
            cnt += 1
print('Number of pathes longer 2: {}'.format(cnt))