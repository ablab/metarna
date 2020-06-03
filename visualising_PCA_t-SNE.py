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
df = pd.concat([X, pca_df, pgm], axis=1)

print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

plt.figure(figsize=(16,10))
pca_plt = sns.scatterplot(
    x="pca_1", y="pca_2",
    hue="regular",
    palette=sns.color_palette("hls", 10),
    data=df,
    legend="full",
    alpha=0.3
)
pca_plt.figure.savefig(os.path.join(sys.argv[3], "pca.png"))