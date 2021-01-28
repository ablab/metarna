from __future__ import print_function

import os, sys

import time

import numpy as np

import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns

from gfa_parser import gfa_to_G, get_one_type_gfa, one_type_gfa_to_df
from spaligner_parser import spaligner_to_df_not_ss


# Coloring using db
# Transcript names define the cluster (i.e. color) of node and all its persons
# Here we don't know how transcripts correspond to persons so can't identify their colors
# because the input graph is regular
def db_coloring(spaligner_ground_truth_tsv, G):
    tsv_df = spaligner_to_df_not_ss(spaligner_ground_truth_tsv, G)
    # Split path column into multiple rows
    new_df = pd.DataFrame(tsv_df['path of the alignment'].str.replace(';', ',').str.split(',').tolist(),
                          index=tsv_df['sequence name']).stack()
    new_df = new_df.reset_index([0, 'sequence name'])
    new_df.columns = ['ground_truth', 'initial_node']
    # Generate set of sequence names for each node with orientation
    db_colors = new_df.groupby('initial_node')['ground_truth'].apply(set).apply(' '.join)
    return db_colors

def persona_coloring(persona_clustering_tsv):
    # Coloring using persona graph clustering
    # Now we know colors for persons separately (not for initial graph nodes)
    # But it isn't ground truth since depends on clustering quality
    persona_colors = pd.Series()
    with open(persona_clustering_tsv, 'r') as fin:
        num_cluster = 0
        for line in fin:
            personas = line.strip().split(',')
            curr = pd.Series([num_cluster] * len(personas), index=personas)
            persona_colors = persona_colors.append(curr, verify_integrity=True)
            num_cluster += 1
    return persona_colors

# Coloring using SPAdes gfa
# Transcript (path) names define the cluster (i.e. color) of node and all its persons
def spades_coloring(gfa, outdir):
    p_gfa = get_one_type_gfa(gfa, 'P', outdir)
    p_gfa_df = one_type_gfa_to_df(p_gfa)
    os.remove(p_gfa)
    colors = pd.DataFrame(p_gfa_df.SegmentNames.str.split(',').tolist(), index=p_gfa_df.PathName).stack()
    colors = colors.reset_index()[[0, 'PathName']]
    colors.columns = ['SegmentNames', 'PathName']
    # To distinguish forward and reverse complement transcript colors
    colors['PathName'] = colors['PathName'].apply(lambda p: "{}+".format(p))
    # Colors for reverse complement nodes
    # since pathes (and links) in gfa includes only one of them (forward or rc)
    rc_colors = colors.applymap(lambda s: s.translate(str.maketrans({'+': '-', '-': '+'})))
    spades_colors = pd.concat([colors, rc_colors], axis=0).set_index('SegmentNames')
    spades_colors = spades_colors.groupby('SegmentNames')['PathName'].apply(set).apply(' '.join)
    return spades_colors

def do_PCA(X):
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(X.values)

    pca_df = pd.DataFrame({'pca_1': pca_result[:, 0],
                           'pca_2': pca_result[:, 1],
                           'pca_3': pca_result[:, 2]},
                          index=X.index)
    print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
    return pca_df

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

def do_umap(X, n_neighbors):
    time_start = time.time()

    umap = UMAP(n_neighbors=n_neighbors, verbose=True)
    umap_result = umap.fit_transform(X.values)

    print('UMAP done! Time elapsed: {} seconds'.format(time.time() - time_start))

    umap_df = pd.DataFrame({'umap_{}_1'.format(n_neighbors): umap_result[:, 0],
                            'umap_{}_2'.format(n_neighbors): umap_result[:, 1]},
                           index=X.index)
    return umap_df

def plot_umap(df, color_col, n_neighbors, outdir):
    plt.figure(figsize=(16, 10))
    umap_plt = sns.scatterplot(
        x="umap_{}_1".format(n_neighbors),
        y="umap_{}_2".format(n_neighbors),
        hue=color_col,
        palette=sns.color_palette("hls", df[color_col].nunique()),
        data=df,
        legend=None,
        # alpha=0.3
    )
    plt.title('n_neighbors = {}'.format(n_neighbors))
    umap_plt.figure.savefig(os.path.join(outdir, "umap.{}.{}.png".format(color_col, n_neighbors)))

# persona_embedding.tsv persona_graph_mapping.tsv node_to_db.tsv persona_clustering.tsv outdir
def visualize_embedding(embedding_df, persona_to_node_tsv, spaligner_ground_truth_tsv, p_clustering_tsv, gfa, G, outdir):
    persona_to_node = pd.read_csv(persona_to_node_tsv, sep=' ',
                                  header=None, index_col=0,
                                  names=['initial_node'])
    df = pd.concat([embedding_df, persona_to_node], axis=1)

    # Coloring using db
    node_colors = db_coloring(spaligner_ground_truth_tsv, G)
    df = df.join(node_colors, on='initial_node')
    # Colorize nodes without pathes in red
    df['ground_truth'] = df['ground_truth'].fillna('0')

    # Coloring using SPAdes pathes
    spades_colors = spades_coloring(gfa, outdir)
    df = df.join(spades_colors, on='initial_node').fillna('0')

    # Coloring using persona graph clustering
    persona_colors = persona_coloring(p_clustering_tsv)
    df = pd.concat([df, persona_colors.to_frame(name='persona_color')], axis=1)

    sns.pairplot(df, vars=embedding_df.keys()).savefig(os.path.join(outdir, "pairplot.png"))

    # PCA
    pca_df = do_PCA(embedding_df)
    df = pd.concat([df, pca_df], axis=1)
    plot_pca_2d(df, 'ground_truth', outdir)
    # plot_pca_3d(df, 'ground_truth')
    plot_pca_2d(df, 'persona_color', outdir)
    # plot_pca_3d(df, 'persona_color')
    plot_pca_2d(df, 'PathName', outdir)

    # T-SNE
    X_subset, df_subset = get_subset(embedding_df, df, 10000)
    # pca_df = do_PCA(X_subset)
    tsne_df = do_t_SNE(X_subset)
    df_subset = pd.concat([df_subset, tsne_df], axis=1)
    plot_t_SNE(df_subset, 'ground_truth', outdir)
    plot_t_SNE(df_subset, 'persona_color', outdir)
    plot_t_SNE(df_subset, 'PathName', outdir)

    # UMAP
    # plot_umap(df_subset, 'ground_truth', 15, outdir)
    # plot_umap(df_subset, 'persona_color', 15, outdir)
    for n in (2, 5, 10, 20, 50, 100, 200):
        umap_df = do_umap(X_subset, n)
        df_subset = pd.concat([df_subset, umap_df], axis=1)
        plot_umap(df_subset, 'PathName', n, outdir)
        plot_umap(df_subset, 'ground_truth', n, outdir)


def main():
    from nxG2clusters import get_total_emb

    outdir = sys.argv[1]
    gfa = sys.argv[2]
    k = int(sys.argv[3])

    persona_to_node_tsv = os.path.join(outdir, 'persona_graph_mapping.tsv')
    node_to_db_tsv = os.path.join(outdir, 'node_to_db.tsv')
    p_clustering_tsv = os.path.join(outdir, 'persona_clustering.tsv')
    p_emb_tsv = os.path.join(outdir, 'persona_embedding.clear.tsv')
    features_tsv = os.path.join(outdir, 'features.tsv')

    tot_emb_df = get_total_emb(p_emb_tsv, features_tsv, persona_to_node_tsv)

    visualize_embedding(tot_emb_df, persona_to_node_tsv, node_to_db_tsv, p_clustering_tsv, gfa, gfa_to_G(gfa, k), outdir)


if __name__ == '__main__':
    main()