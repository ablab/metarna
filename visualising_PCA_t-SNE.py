# visualising_PCA_t-SNE.py persona_embedding.tsv persona_graph_mapping.tsv outdir

from __future__ import print_function

import os, sys

import time

import numpy as np

import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns

X = pd.read_csv(sys.argv[1], sep=' ', header=None, index_col=0, skiprows=1, skipfooter=10)
pgm = pd.read_csv(sys.argv[2], sep=' ', header=None, index_col=0, names=['regular'])

pca = PCA(n_components=3)
pca_result = pca.fit_transform(X.values)

pca_df = pd.DataFrame({'pca_1': pca_result[:, 0],
                       'pca_2': pca_result[:, 1],
                       'pca_3': pca_result[:, 2]}, index=X.index)
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

df = pd.concat([X, pca_df, pgm], axis=1)

# if clusters overlap then 0 else cluster number
# path1 930004-,278546-,36185+,278990+,283130+,352975-,37703+
# path2 930004-,239212-,36185+,365256-,283130+,352975-,37703+
clusters_map = {'930004-': '0', '36185+': '0', '283130+': '0', '352975-': '0', '37703+': '0',
                '278546-': '1', '278990+': '1',
                '239212-': '2', '365256-': '2',
                '2326645-': '-1'}
df['cluster'] = df['regular'].map(clusters_map)

plt.figure(figsize=(16,10))
pca_plt = sns.scatterplot(
    x="pca_1", y="pca_2",
    hue="cluster",
    palette=sns.color_palette("hls", 4),
    data=df,
    legend="full",
    # alpha=0.3
)
pca_plt.figure.savefig(os.path.join(sys.argv[3], "pca.png"))