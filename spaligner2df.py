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