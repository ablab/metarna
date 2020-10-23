import community

# from {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 1 ... }
# to [[0, 1, 2, 3], [4, 5, 6], ... ]
def clusters_dict_to_list(c_dict):
    c_list = []
    for cluster in set(c_dict.values()):
        c_list.append([nodes for nodes in c_dict.keys() if c_dict[nodes] == cluster])
    return c_list

def best_partition_cov_diff(G, weight='cov_diff'):
    c_dict = community.best_partition(G, weight=weight)
    c_list = clusters_dict_to_list(c_dict)
    return c_list

def best_partition_long_reads(G, weight='num_long_reads'):
    c_dict = community.best_partition(G, weight=weight)
    c_list = clusters_dict_to_list(c_dict)
    return c_list