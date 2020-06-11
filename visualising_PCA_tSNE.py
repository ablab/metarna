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


def do_PCA(X):
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(X.values)

    pca_df = pd.DataFrame({'pca_1': pca_result[:, 0],
                           'pca_2': pca_result[:, 1],
                           'pca_3': pca_result[:, 2]},
                          index=X.index)
    print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
    return pca_df

def persona_coloring(persona_clustering_tsv):
    # Coloring using persona graph clustering
    # Now we know colors for persons separately (not for initial graph nodes)
    # But it isn't ground truth since depends on clustering quality
    persona_colors = pd.Series()
    with open(persona_clustering_tsv, 'r') as fin:
        num_cluster = 0
        for line in fin:
            personas = [int(p) for p in line.split()]
            curr = pd.Series([num_cluster] * len(personas), index=personas)
            persona_colors = persona_colors.append(curr, verify_integrity=True)
            num_cluster += 1
    return persona_colors

# PCA
def plot_pca_2d(df, color_col, outdir):
    plt.figure(figsize=(16, 10))
    pca_plt = sns.scatterplot(
        x="pca_1", y="pca_2",
        hue=color_col,
        palette=sns.color_palette("hls", df[color_col].nunique()),
        data=df,
        legend=None,
        # alpha=0.3
    )
    pca_plt.figure.savefig(os.path.join(outdir, "pca_2d.{}.png".format(color_col)))

def plot_pca_3d(df, color_col, outdir):
    ax = plt.figure(figsize=(16, 10)).gca(projection='3d')
    ax.scatter(
        xs=df["pca_1"],
        ys=df["pca_2"],
        zs=df["pca_3"],
        c=df[color_col],
        cmap='tab10'
    )
    ax.set_xlabel('pca-one')
    ax.set_ylabel('pca-two')
    ax.set_zlabel('pca-three')
    plt.savefig(os.path.join(outdir, "pca_3d.{}.png".format(color_col)))

def get_subset(X, df, N=10000):
    X_subset = X.sample(min(df.shape[0], N))
    df_subset = df.loc[X_subset.index, :]

    return X_subset, df_subset

# Since t-SNE scales quadratically in the number of objects N,
# its applicability is limited to data sets with only a few thousand input objects.
def do_t_SNE(X):
    time_start = time.time()

    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_result = tsne.fit_transform(X.values)

    print('t-SNE done! Time elapsed: {} seconds'.format(time.time() - time_start))

    tsne_df = pd.DataFrame({'tsne_1': tsne_result[:, 0],
                            'tsne_2': tsne_result[:, 1]},
                           index=X.index)
    return tsne_df

def plot_t_SNE(df, color_col, outdir):
    plt.figure(figsize=(16, 10))
    t_SNE_plt = sns.scatterplot(
        x="tsne_1", y="tsne_2",
        hue=color_col,
        palette=sns.color_palette("hls", df[color_col].nunique()),
        data=df,
        legend=None,
        # alpha=0.3
    )
    t_SNE_plt.figure.savefig(os.path.join(outdir, "t-SNE.{}.png".format(color_col)))

# persona_embedding.tsv persona_graph_mapping.tsv node_to_db.tsv persona_clustering.tsv outdir
def visualize_embedding(p_emb_tsv, persona_to_node_tsv, node_to_db_tsv, p_clustering_tsv, outdir):
    p_emb = pd.read_csv(p_emb_tsv, sep=' ', header=None, index_col=0, skiprows=1)

    persona_to_node = pd.read_csv(persona_to_node_tsv, sep=' ', header=None, index_col=0, names=['initial_node'])

    pca_df = do_PCA(p_emb)

    df = pd.concat([p_emb, pca_df, persona_to_node], axis=1)

    # Coloring using db
    # Transcript names define the cluster (i.e. color) of node and all its persons
    # Here we don't know how transcripts correspond to persons so can't identify their colors
    # because the input graph is regular
    node_colors = pd.read_csv(node_to_db_tsv, sep='\t', header=None, index_col=0, names=['ground_truth'])
    df = df.join(node_colors, on='initial_node')
    # colorize nodes without pathes in red
    df['ground_truth'] = df['ground_truth'].fillna('0')

    persona_colors = persona_coloring(p_clustering_tsv)
    df = pd.concat([df, persona_colors.to_frame(name='persona_color')], axis=1)

    plot_pca_2d(df, 'ground_truth', outdir)
    # plot_pca_3d(df, 'ground_truth')
    plot_pca_2d(df, 'persona_color', outdir)
    # plot_pca_3d(df, 'persona_color')

    X_subset, df_subset = get_subset(p_emb, df, 10000)
    # pca_df = do_PCA(X_subset)
    tsne_df = do_t_SNE(X_subset)
    df_subset = pd.concat([df_subset, tsne_df], axis=1)

    plot_t_SNE(df_subset, 'ground_truth', outdir)
    plot_t_SNE(df_subset, 'persona_color', outdir)