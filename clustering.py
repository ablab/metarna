import community

# from {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 1 ... }
# to [[0, 1, 2, 3], [4, 5, 6], ... ]
def clusters_dict_to_list(c_dict):
    c_list = []
    for cluster in set(c_dict.values()):
        c_list.append([nodes for nodes in c_dict.keys() if c_dict[nodes] == cluster])
    return c_list

def clusters_list_to_dict(c_list):
    c_dict = {}
    for i, com in enumerate(c_list):
        for node in com:
            c_dict[node] = i
    return c_dict

def best_partition(G, weight='weight'):
    c_dict = community.best_partition(G, weight=weight)
    c_list = clusters_dict_to_list(c_dict)
    return c_list