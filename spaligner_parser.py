import pandas as pd

from Bio.Seq import Seq


def spaligner_to_df(tsv):
    tsv_df = pd.read_csv(tsv, sep="\t", names=['sequence name',
                                               'start position of alignment on sequence',
                                               'end position of alignment on sequence',
                                               'start position of alignment on the first edge of the Path',
                                               'end position of alignment on the last edge of the Path',
                                               'sequence length',
                                               'path of the alignment',
                                               'lengths of the alignment on each edge of the Path respectively',
                                               'sequence of alignment Path'])
    return tsv_df

def rc_smth_of_the_alignment(path_str):
    pathes = path_str.split(';')
    rc_pathes = []
    for path in pathes:
        compl_path = path.translate(str.maketrans({'+': '-', '-': '+'}))
        rc_path = ','.join(reversed(compl_path.split(',')))
        rc_pathes.append(rc_path)
    rc_path_str = ';'.join(reversed(rc_pathes))
    return rc_path_str

def spaligner_to_df_not_ss(tsv, G):
    tsv_df = spaligner_to_df(tsv)
    tsv_df_rc = tsv_df.copy()
    tsv_df_rc['sequence name'] = tsv_df['sequence name'].astype(str) + '_rc'

    path = tsv_df['path of the alignment'].str.replace(';', ',').str.split(',')
    start_node = path.str[0]
    end_node = path.str[-1]
    s_pos = tsv_df['start position of alignment on sequence'].astype(str).str.split(',').str[0].astype(int)
    e_pos = tsv_df['end position of alignment on sequence'].astype(str).str.split(',').str[-1].astype(int)
    s_len = start_node.apply(lambda x: G.nodes[x]['len']).astype(int)
    e_len = end_node.apply(lambda x: G.nodes[x]['len']).astype(int)
    tsv_df_rc['start position of alignment on the first edge of the Path'] = e_len - e_pos - G.graph['k']
    tsv_df_rc['end position of alignment on the last edge of the Path'] = s_len - s_pos - G.graph['k']

    tsv_df_rc['path of the alignment'] = \
        tsv_df['path of the alignment'].apply(rc_smth_of_the_alignment)
    tsv_df_rc['lengths of the alignment on each edge of the Path respectively'] = \
        tsv_df['lengths of the alignment on each edge of the Path respectively'].apply(rc_smth_of_the_alignment)
    tsv_df_rc['sequence of alignment Path'] = \
        tsv_df['sequence of alignment Path'].map(lambda s: Seq(s).reverse_complement())
    return pd.concat([tsv_df, tsv_df_rc])

def spaligner_to_clustering_tsv(spaligner_tsv, clustering_tsv, G, min_clusters_size=2):
    tsv_df = spaligner_to_df_not_ss(spaligner_tsv, G)
    # tsv_df['path of the alignment'] = tsv_df['path of the alignment'].str.replace(',', ' ')
    pathes = pd.Series(tsv_df['path of the alignment']).str.replace(';', ',')
    filtered_pathes = pathes[pathes.str.split(',').str.len() >= min_clusters_size]
    clusters = filtered_pathes.drop_duplicates()
    clusters.to_csv(clustering_tsv, sep='\t', header=False, index=False)
    return clustering_tsv