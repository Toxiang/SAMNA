import networkx as nx
import random

class G:
    max_weight = -1.0
    Nk = 0


# kset = []


def find_clique(g, k, included):
    # kset.clear()
    G.Nk = k
    node_best_weights = {}
    max_clique = []
    add_for_k = {}
    G.max_weight = -1.0
    all_nodes = []
    clique = []
    kset = []

    not_hide_nodes = []

    for v1 in included:
        clique.append(v1)
        max_clique.append(v1)
        not_hide_nodes.append(v1)
        if v1[0:2] not in kset:
            kset.append(v1[0:2])
    for v2 in g.nodes():         # only left node neighbour to clique
        if not clique.__contains__(v2):
            add = True
            for i in clique:
                res = False
                if g.has_edge(v2, i):
                    res = True
                if res is False:
                    add = False
            if add:
                all_nodes.append(v2)
                not_hide_nodes.append(v2)
    g1 = nx.Graph(g.subgraph(not_hide_nodes))    # graph of hiden node
    l = []
    for v3 in g1.nodes():
        dic = {}
        for n in g1.neighbors(v3):
            if g1.has_edge(v3, n):
                weight = g1.get_edge_data(v3, n)['weight']
            net = n[0:2]
            if not dic.__contains__(net):
                dic[net] = weight
            elif dic[net] < weight:
                dic[net] = weight
            l.append(0.0 - weight)
        node_best_weights[v3] = dic.copy()
    add_for_k[1] = 0.0
    for i2 in range(2, k):
        add_for_k[i2] = add_for_k[i2-1]
        for j in range(1, i2):
            if len(l) > 1:
                add_for_k[i2] -= min(l)
                l.remove(min(l))
                l.remove(min(l))
    if len(all_nodes) > 0:
        g1, all_nodes, max_clique = scan(0, len(kset)+1, g1, clique, 0.0, all_nodes, node_best_weights, add_for_k,
                                         max_clique)

    #if len(max_clique) < 3:

    result = []
    if len(max_clique) > 1:
        for v4 in g:
            if v4 in max_clique:
                result.append(v4)
    return result


def scan(l_int, d, g1, clique, clique_weight, neighbours, node_best_weights, add_for_k, max_clique):
    while len(neighbours) > 0:
        v_scan = neighbours.pop(0)
        intersection = []
        kset_scan = []
        new_clique = clique.copy()
        new_clique_weight = clique_weight
        cont = True
        updated = False
        add_new_node = False
        for i_1 in clique:
            if g1.has_edge(i_1, v_scan):
                new_clique_weight += g1.get_edge_data(i_1, v_scan)['weight']
                add_new_node = True
            else:
                add_new_node = False
                # print(f"Alert  {i_1}  {v_scan} ")
        new_clique.append(v_scan)   # add node even if can not be a clique?
        if new_clique_weight >= G.max_weight:
            max_clique = new_clique.copy()
            G.max_weight = new_clique_weight
            updated = True
        if d < G.Nk:
            for n in neighbours:
                w = -1.0
                if g1.has_edge(n, v_scan):
                    w = g1.get_edge_data(n, v_scan)['weight']
                if w > 0.0:
                    intersection.append(n)
                    if n[0:2] not in kset_scan:
                        kset_scan.append(n[0:2])
            intersection_k = len(kset_scan)
            if l_int == 0:
                g1.remove_node(v_scan)
            if intersection_k == 0:
                cont = False
            elif updated is False:
                to_add = add_for_k[intersection_k]
                for i_scan in new_clique:
                    for net in kset_scan:
                        # print(i)
                        # if node_best_weights[i_scan].__contains__(net):
                        to_add += node_best_weights[i_scan][net]
                if (new_clique_weight + to_add) < G.max_weight:
                    cont = False
            if cont:
                g1, neighbours, max_clique = scan(l_int+1, d+1, g1, new_clique, new_clique_weight,
                                                  intersection, node_best_weights, add_for_k, max_clique)
    return g1, neighbours, max_clique


def take_second(elem):
    return elem[1]
