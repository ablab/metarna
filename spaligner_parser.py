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

# Lines such as the following
# NODE_746_length_903_cov_38.166276_g539_i1       0       903     355     80      903
# 42629-,22255+,40519-,38909+,41393+,38139+       0,285,90,54,394,80      AGGGTGTGGCAGAGGCAG...
num_zero_alignment_lines = 0
num_zero_alignments = 0

def remove_zero_alignments(series_obj):
    global num_zero_alignment_lines
    global num_zero_alignments

    # TODO Need to recalculate start and end after zero alignments removing
    # TODO since there are ; and else cases
    # series_obj['start position of alignment on sequence'] = '*'
    # series_obj['end position of alignment on sequence'] = '*'
    # series_obj['start position of alignment on the first edge of the Path'] = '*'
    # series_obj['end position of alignment on the last edge of the Path'] = '*'

    pathes = series_obj['path of the alignment'].split(';')
    lengths = series_obj['lengths of the alignment on each edge of the Path respectively'].split(';')
    new_pathes = []
    new_lengths = []
    for i_path, path in enumerate(pathes):
        nodes = path.split(',')
        lengths_for_nodes = lengths[i_path].split(',')
        nonzero_ind = [i for i, l in enumerate(lengths_for_nodes) if l != 0]
        n_diff = len(lengths_for_nodes) - len(nonzero_ind)
        num_zero_alignments += n_diff
        if n_diff > 0:
            num_zero_alignment_lines += 1
        new_pathes.append(','.join([nodes[i] for i in nonzero_ind]))
        new_lengths.append(','.join([lengths_for_nodes[i] for i in nonzero_ind]))
    series_obj['path of the alignment'] = ';'.join(new_pathes)
    series_obj['lengths of the alignment on each edge of the Path respectively'] = ';'.join(new_lengths)
    return series_obj

def remove_zero_length_alignments(tsv_in, tsv_out):
    global num_zero_alignment_lines
    global num_zero_alignments

    tsv_df = spaligner_to_df(tsv_in)
    filtered_df = tsv_df.apply(remove_zero_alignments, axis=1)
    filtered_df.to_csv(tsv_out, sep='\t', header=False, index=False)
    print(tsv_in)
    print('Number of lines containing zero length alignments: {}'.format(num_zero_alignment_lines))
    print('Number of zero length alignments: {}\n'.format(num_zero_alignments))
    return tsv_out


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

    start_node = tsv_df['path of the alignment'].str.replace(';', ',').str.split(',').str[0]
    end_node = tsv_df['path of the alignment'].str.replace(';', ',').str.split(',').str[-1]
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
    tsv_df['path of the alignment'] = tsv_df['path of the alignment'].str.replace(',', ' ')
    tsv_df['path of the alignment'] = tsv_df['path of the alignment'].str.replace(';', ' ')
    tsv_df = tsv_df[tsv_df['path of the alignment'].str.split(' ').str.len() >= min_clusters_size]
    tsv_df.to_csv(clustering_tsv,
                  columns=['path of the alignment'],
                  sep='\t', header=False, index=False)
    return clustering_tsv