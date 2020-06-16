import pandas as pd


def spaligner_to_df(tsv):
    tsv_df = pd.read_csv(tsv, sep="\t", names=['sequence name',
                                               'start position of alignment on sequence',
                                               'end position of  alignment on sequence',
                                               'start position of alignment on the first edge of the Path',
                                               'end position of alignment on the last edge of the Path',
                                               'sequence length',
                                               'path of the alignment',
                                               'lengths of the alignment on each edge of the Path respectively',
                                               'sequence of alignment Path'])
    return tsv_df

def spaligner_to_clustering_tsv(spaligner_tsv, clustering_tsv):
    tsv_df = spaligner_to_df(spaligner_tsv)
    tsv_df['path of the alignment'] = tsv_df['path of the alignment'].str.replace(',', ' ')
    tsv_df['path of the alignment'] = tsv_df['path of the alignment'].str.replace(';', ' ')
    tsv_df.to_csv(clustering_tsv,
                  columns=['path of the alignment'],
                  sep='\t', header=False, index=False)
    return clustering_tsv