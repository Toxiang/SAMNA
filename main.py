import networkx as nx
import sys
import numpy as np
import time
import math
import datetime
import random
import simulatedAnnealing as SA
import maxclique as mc


class C:
    k = 0

def readnetwork(datafile):
    f = open(datafile, 'r')
    lines = f.readlines()
    print("datafile opened.")
    #global k
    C.k = int(lines[0])  # number of networks
    # print(k)
    print("reading networks file......")
    for i in range(1, C.k + 1):
        netfile = lines[i].strip()
        print("opening  {netfile}  ...".format(netfile=netfile))
        try:
            with open(netfile, 'r') as f1:
                nn = 0
                ne = 0
                for j in f1:
                    lineArr = j.strip().split()
                    node1 = lineArr[0]
                    node2 = lineArr[1]
                    if ppi.has_node(node1):
                        node1 = lineArr[0]
                    else:
                        ppi.add_node(node1)
                        sim_ppi.add_node(node1)
                        node_prio.append(node1)
                        nn = nn + 1
                    if ppi.has_node(node2):
                        node2 = lineArr[1]
                    else:
                        ppi.add_node(node2)
                        sim_ppi.add_node(node2)
                        node_prio.append(node2)
                        nn = nn + 1
                    if not ppi.has_edge(node1, node2):
                        ppi.add_edge(node1, node2)
                        ne = ne + 1
                        # node_pair = [node1, node2]
                        # EDGE_weights[node_pair] = 1
                # print(ppi.node())
                print("There are {nn} nodes and {ne} edge.".format(nn=nn,ne=ne))
        except FileNotFoundError:
            print("can not find {netfile} ppi file!\n".format(netfile=netfile))
    best_weight = {}
    node_Best = {}
    best_sim = {}
    no_self_sim_list = []

    print("reading BLAST similiarity file......")
    for i1 in lines[C.k + 1:]:
        netfile = i1.strip()
        print("opening {netfile}. ".format(netfile=netfile))
        try:
            with open(netfile, 'r') as f2:
                for j1 in f2:
                    lineArr = j1.strip().split()
                    sim_node1 = lineArr[0]
                    sim_node2 = lineArr[1]
                    w = float(lineArr[2])
                    # print (sim_node1[0:2])
                    if ppi.has_node(sim_node1):
                        if ppi.has_node(sim_node2):
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)
                        else:
                            ppi.add_node(sim_node2)
                            sim_ppi.add_node(sim_node2)
                            node_prio.append(sim_node2)
                            nn = nn + 1
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)

                    else:
                        ppi.add_node(sim_node1)
                        sim_ppi.add_node(sim_node1)
                        node_prio.append(sim_node1)
                        nn = nn + 1
                        if ppi.has_node(sim_node2):
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)

                        else:
                            ppi.add_node(sim_node2)
                            sim_ppi.add_node(sim_node2)
                            node_prio.append(sim_node2)
                            nn = nn + 1
                            if sim_node1 != sim_node2:
                                sim_ppi.add_edge(sim_node1, sim_node2, weight=w)


                    # visit and find each node's max weight of connected node
                    if best_weight.__contains__(sim_node1):
                        if w > best_weight[sim_node1]:
                            best_weight[sim_node1] = w
                    else:
                        best_weight[sim_node1] = w
                    if best_weight.__contains__(sim_node2):
                        if w > best_weight[sim_node2]:
                            best_weight[sim_node2] = w
                    else:
                        best_weight[sim_node2] = w

                    if not best_sim.__contains__(sim_node1):
                        best_sim[sim_node1] = w
                    elif w > best_sim[sim_node1]:
                        best_sim[sim_node1] = w
                    if not best_sim.__contains__(sim_node2):
                        best_sim[sim_node2] = w
                    elif w > best_sim[sim_node2]:
                        best_sim[sim_node2] = w
                    if sim_node1[0:2] != sim_node2[0:2]:
                        no_self_sim_list.append((sim_node1, sim_node2))
                        nodepair = (sim_node1, sim_node2)
                        simmap[nodepair] = w

                    d1 = {}
                    d2 = {}
                    if node_Best.__contains__(sim_node1):
                        d1 = node_Best[sim_node1]
                    if node_Best.__contains__(sim_node2):
                        d2 = node_Best[sim_node2]
                    if not d1.__contains__(sim_node2[0:2]):
                        d1[sim_node2[0:2]] = 0.0
                    if not d2.__contains__(sim_node1[0:2]):
                        d2[sim_node1[0:2]] = 0.0

                    if d1[sim_node2[0:2]] < w:
                        d1[sim_node2[0:2]] = w
                    if d2[sim_node1[0:2]] < w:
                        d2[sim_node1[0:2]] = w

                    node_Best.setdefault(sim_node1, d1)
                    node_Best.setdefault(sim_node2, d2)
                f2.close()
        except FileNotFoundError:
            print("can not find {netfile} BLAST similarity file!\n".format(netfile=netfile))
    f.close()
    print("Finished reading BLAST similarities.")
    for tt in no_self_sim_list:
        simmap[tt] = simmap[tt] / (pow(best_sim[tt[0]], 0.5) * pow(best_sim[tt[1]], 0.5))
    no_self_sim_list.clear()
    best_sim.clear()
    # print (best_weight)
    print("In total, There are {nn} nodes and {ne} edges.".format(nn=sim_ppi.number_of_nodes(),ne=sim_ppi.number_of_edges()))
    # filter
    print("Filtering BLAST scores.")
    # print(sim_ppi.edges(data='weight'))
    remove_node = []


    for node1, node2, fil_w in sim_ppi.edges(data='weight'):
        # for node1, nbrs in sim_ppi.adjacency():
        #     for node2, e_data in nbrs.items():
        #         fil_w = e_data['weight']
        if fil_w < beta * node_Best[node1][node2[0:2]] or fil_w < beta * node_Best[node2][node1[0:2]]:
            remove_node.append((node1, node2))
    for node_pair in remove_node:
        sim_ppi.remove_edge(node_pair[0], node_pair[1])
    weight_list = []
    for node1, node2, w in sim_ppi.edges(data='weight'):
        weight_list.append(w)
    x = np.array(weight_list)
    max_w = x.max()
    min_w = x.min()
    ranges = max_w - min_w
    print("There are {nn} ppi nodes and {ne} ppi edges.".format(nn=ppi.number_of_nodes(),ne=ppi.number_of_edges()))
    print(
        "After all, There are {nn} sim_ppi nodes and {ne} sim_ppi edges.".format(nn=sim_ppi.number_of_nodes(), ne=sim_ppi.number_of_edges()))


def construct_star(node):
    neighbors = list(sim_ppi.neighbors(node))
    total = 0.0
    matchneighbor = [node]
    for nei in neighbors:
        if not node_to_cluster_map.__contains__(nei):
            matchneighbor.append(nei)
    star = nx.Graph(sim_ppi.subgraph(matchneighbor))
    return star

def createclique():
    for prio in sim_ppi.nodes():
        print(prio)
        prio_graph = construct_star(prio)
        nodes = []
        nodes.append(prio)
        prio_clique = mc.find_clique(prio_graph, C.k, nodes)
        if len(prio_clique)>1:
            NODEs_clique[prio] = prio_clique
            NODEs_Weight[prio] = CW(prio_clique)
            candidate_clique.append(prio_clique)

def CW(clique):
    result = 0.0
    for u in clique:
        for v in clique:
            if sim_ppi.has_edge(u,v):
                result = result + sim_ppi.get_edge_data(u,v)['weight']
    return result/2

def simuldated_annealing(candidate):
    t_Max = 100  #initiate temperature
    t_Min = 10   #minimum value of terperature
    k_Max = 100
    k = 100   #number of iterations
    s = 0.005
    n_Max = 2000

    k = 0
    t = t_Max
    step = (t_Max - t_Min) / k_Max
    while k <= k_Max and len(candidate)>0:
        if len(candidate)==0:
            break
        t = t - step
        beta_s = 1.0 / (s * t)
        for i in range(0, n_Max):
            if len(candidate)==0:
                break
            delta = 0.0
            print(len(candidate),k)
            select = random.choice(candidate)
            overlapValue, delta, clu = isOverlapped_n(select)
            if not overlapValue:
                add(select, clu)
                candidate.remove(select)
            elif overlapValue == 3:
                candidate.remove(select)
                new_cluster = []
                for u in select:
                    if not node_to_cluster_map.__contains__(u) and not labelnode.__contains__(u):
                        new_cluster.append(u)
                add(new_cluster, -1)
                # candidate.append(new_cluster)
            elif overlapValue == 4:
                candidate.remove(select)
            elif overlapValue ==5:
                candidate.remove(select)
            elif overlapValue ==6:
                candidate.remove(select)
                new_cluster=[]
                for u in select:
                    if not node_to_cluster_map.__contains__(u) and not labelnode.__contains__(u):
                        new_cluster.append(u)
                candidate.append(new_cluster)
            else:
                if delta > 0 or random.random() < math.exp(beta_s*delta):
                    if(updateState(select, delta, overlapValue, clu)):
                        candidate.remove(select)
        k=k+1

def isOverlapped_n(select):
    overlapValue = 0
    length = len(select)
    deltaScore = [0] * 10
    clu_index = [-1] * 10
    align_node = {}
    num = 0
    node_cluster_map = {}
    for node in select:
        if node in labelnode:
            num=num+1
            clusters_i = node_to_cluster_map[node]
            if not node_cluster_map.__contains__(clusters_i):
                node_cluster_map[clusters_i] = 1
            else:
                node_cluster_map[clusters_i] = node_cluster_map[clusters_i]+1
    if node_cluster_map:
        max_num = max(node_cluster_map.values())
        if max_num > 0.5*length and max_num<length:
            print("here right!")
            cl_index = max(node_cluster_map, key = node_cluster_map.get)
            # clu_index[overlapValue] = cl_index
            com_cluster = clusters[cl_index].copy()
            new_cluster = com_cluster.copy()
            net = []
            for u in com_cluster:
                if not net.__contains__(u[0:2]):
                    net.append(u[0:2])
            if isCombine_n(com_cluster,select):
                for u in select:
                    if not net.__contains__(u[0:2]):
                        new_cluster.append(u)
                overlapValue = 1
                deltaScore[1] = getDeltaScore(new_cluster, com_cluster, cl_index)
                clu_index[1] = cl_index
            else:
                for u in select:
                    if not com_cluster.__contains__(u) and net.__contains__(u[0:2]):
                        for v in com_cluster:
                            if v[0:2]==u[0:2]:
                                new_cluster.remove(v)
                                new_cluster.append(u)
                                break
                overlapValue = 2
                deltaScore[2] = getDeltaScore(new_cluster, com_cluster, cl_index)
                clu_index[2] = cl_index
        elif num <= 0.25*length:
            new_cluster = []
            for clu in node_cluster_map:
                node_num = node_cluster_map[clu]
                com_cluster = clusters[clu]
                for u in select:
                    if not com_cluster.__contains__(u):
                        new_cluster.append(u)
            overlapValue = 3
            empty = []
            deltaScore[3] = getDeltaScore(new_cluster, empty, 0)
        elif max_num == length:
            overlapValue = 4
        elif num == length and max_num < length:
            overlapValue = 5
        else:
            overlapValue = 6
    myDeltaScore = deltaScore[overlapValue]
    cluster_index = clu_index[overlapValue]
    return overlapValue, myDeltaScore, cluster_index

def isOverlapped(select):
    overlapValue = 0
    protein1 = select[0]
    protein2 = select[1]
    species1 = protein1[0:2]
    species2 = protein2[0:2]
    deltaScore = [0] * 10
    clu_index = [-1] * 10
    index_i = 0
    index_j = 0
    if protein1 in labelnode and protein2 in labelnode:
        c_n1 = cluster_node_map[protein1]
        c_n2 = cluster_node_map[protein2]
        cluster1 = clusters[c_n1]
        cluster2 = clusters[c_n2]
        newmp1 = cluster1.copy()
        newmp2 = cluster2.copy()
        if isCombine(cluster1, protein2, species2):
            newmp1.append(protein2)
            for node in newmp2:
                if node == protein2:
                    newmp2.remove(node)
            deltaScore[1] = getDeltaScore(newmp1, cluster1, c_n1) + getDeltaScore(newmp2, cluster2, c_n2)
            clu_index[1] = c_n1
            index_i = 1
        else:
            for u in newmp1:
                if u[0:2]==species2:
                    newmp1 = [protein2 if node == u else node for node in newmp1]
                    # newmp1 = newmp1.replace(u, protein2)
                    for v in newmp2:
                        if v == protein2:
                            newmp2.remove(v)
            deltaScore[2] = getDeltaScore(newmp1, cluster1, c_n1) + getDeltaScore(newmp2, cluster2, c_n2)
            clu_index[2] = c_n1
            index_i = 2
        newmp1 = cluster1
        newmp2 = cluster2
        if isCombine(cluster2, protein1,species1):
            newmp2.append(protein1)
            for node in newmp1:
                if node == protein1:
                    newmp1.remove(node)
            deltaScore[3] = getDeltaScore(newmp2, cluster2, c_n2) + getDeltaScore(newmp1, cluster1, c_n1)
            clu_index[3] = c_n2
            index_j = 3
        else:
            for u in newmp2:
                if u[0:2] == species1:
                    newmp2 = [protein1 if node == u else node for node in newmp2]
                    #newmp2 = newmp2.replace(u, protein1)
                    for v in newmp1:
                        if v == protein2:
                            newmp1.remove(v)
            deltaScore[4] = getDeltaScore(newmp1,cluster1,c_n1) + getDeltaScore(newmp2, cluster2,c_n2)
            clu_index[4] = c_n2
            index_j = 4
        maxDelta = deltaScore[index_i]
        overlapValue = index_i
        if maxDelta < deltaScore[index_j]:
            maxDelta = deltaScore[index_j]
            overlapValue = index_j
    elif protein1 in labelnode:
        c_n1 = cluster_node_map[protein1]
        cluster1 = clusters[c_n1]
        newmp1 = cluster1.copy()
        if isCombine(cluster1, protein2,species2):
            newmp1.append(protein2)
            deltaScore[5] = getDeltaScore(newmp1, cluster1, c_n1)
            clu_index[5] = c_n1
            overlapValue = 5
        else:
            for u in newmp1:
                if u[0:2] == species2:
                    newmp1 = [protein2 if node == u else node for node in newmp1]
                    #newmp1 = newmp1.replace(u, protein2)
                    deltaScore[6] = getDeltaScore(newmp1, cluster1 ,c_n1)
                    clu_index[6] = c_n1
                    overlapValue = 6
    elif protein2 in labelnode:
        c_n2 = cluster_node_map[protein2]
        cluster2 = clusters[c_n2]
        newmp2 = cluster2.copy()
        if isCombine(cluster2,protein1,species1):
            newmp2.append(protein1)
            deltaScore[7] = getDeltaScore(newmp2, cluster2, c_n2)
            clu_index[7] = c_n2
            overlapValue = 7
        else:
            for u in newmp2:
                if u[0:2] == species1:
                    newmp2 = [protein1 if node == u else node for node in newmp2]
                    # newmp2 = newmp2.replace(u, protein1)
                    deltaScore[8] = getDeltaScore(newmp2, cluster2, c_n2)
                    clu_index[8] = c_n2
                    overlapValue = 8
    myDeltaScore = deltaScore[overlapValue]
    cluster_index = clu_index[overlapValue]
    return overlapValue, myDeltaScore, cluster_index


def isCombine(matchset, protein, species):
    for u in matchset:
        if u[0:2] == species:
            if not sim_ppi.has_edge(u,protein):
                return False
        else:
            return True
    return True

def isCombine_n(matchset1, matchset2):
    net = []
    nodes = []
    for node in matchset1:
        if not net.__contains__(node[0:2]):
            net.append(node[0:2])
    for u in matchset2:
        if not net.__contains__(u[0:2]):
            return True
        elif not matchset1.__contains__(u):
            nodes.append(u)
    for u in nodes:
        flag1 = True
        for v in matchset1:
            if not sim_ppi.has_edge(u,v):
                flag1 = False
        if not flag1:
            return False
    return True


def add(candidate, i):
    if i<0:
        i = len(clusters)
    if len(candidate)>1:
        for node in candidate:
            if not cluster_node_map.__contains__(node) and not labelnode.__contains__(node):
                labelnode.append(node)
                cluster_node_map[node] = i
                node_to_cluster_map[node] = i
        clusters[i] = candidate

def getDeltaScore(newcluster,oldcluster,i):
    score1 = 0.0
    score2 = 0.0
    icq_new, kk_n = calculate_icq(newcluster)
    ciq_new, int_cons = calculate_ciq(newcluster)
    if ciq_new[1] > 0.0:
        score1 = ((1.0 - alpha) * icq_new + alpha * (ciq_new[0] /ciq_new[1]))
    else:
        score1 = (1.0 - alpha) * icq_new
    icq_old, kk_o = calculate_icq(oldcluster)
    ciq_old, int_cons = calculate_ciq(oldcluster)
    if ciq_old[1] > 0.0:
        score2 = ((1.0 - alpha) * icq_old + alpha * (ciq_old[0] / ciq_old[1]))
    else:
        score2 = (1.0 - alpha) * icq_old
    delta = score1 - score2
    return delta

def updateState(select,delta,overlapValue,clu_index):
    old_cluster = clusters[clu_index]
    dic1 = select.copy()
    dic2 = old_cluster.copy()
    net = []
    for u in old_cluster:
        if not net.__contains__(u[0:2]):
            net.append(u[0:2])
    if overlapValue == 1:
        for u in select:
            if not net.__contains__(u[0:2]):
                old_cluster.append(u)
                node_to_cluster_map[u]=clu_index
                labelnode.append(u)
        clusters[clu_index] = old_cluster
    elif overlapValue == 2:
        for u in dic1:
            if not dic2.__contains__(u) and net.__contains__(u[0:2]):
                for v in dic2:
                    if u[0:2]==v[0:2] and labelnode.__contains__(v):
                        old_cluster.remove(v)
                        # del node_to_cluster_map[v]
                        labelnode.remove(v)
                        old_cluster.append(u)
                        node_to_cluster_map[u] = clu_index
                        labelnode.append(u)
                        break
        clusters[clu_index] = old_cluster
    elif overlapValue == 3:
        new_cluster = []
        clu = len(clusters)
        for u in dic1:
            if not node_to_cluster_map.__contains__(u):
                new_cluster.append(u)
                labelnode.append(u)
                node_to_cluster_map[u] = clu
        clusters[clu] = new_cluster
    return True



def calculate_icq(cluster):
    result = 0.0
    sn = 0
    kset_icq = []
    for v in cluster:
        net_name = v[0:2]
        if not kset_icq.__contains__(net_name):
            kset_icq.append(net_name)
        cluster2 = cluster.copy()
        for n in cluster2:
            if v[0:2] != n[0:2]:
                sn += 1
                nodepair1 = (n, v)
                nodepair2 = (v, n)
                if simmap.__contains__(nodepair1):
                    result += simmap[nodepair1]
                elif simmap.__contains__(nodepair2):
                    result += simmap[nodepair2]
                # if sim_ppi.has_edge(v, n):
                #     result = result + sim_ppi.get_edge_data(n, v)['weight']
    if sn != 0:
        result /= sn
    else:
        result = 0.0
    kk = len(kset_icq)
    return result, kk

def calculate_ciq(cluster):
    result = (0.0, 0)
    related_clusters = []
    internal_conserved = 0
    for v in cluster:
        for n in cluster:
            if ppi.has_edge(v, n) and v[0:2] == n[0:2]:
                internal_conserved += 1
        for n in list(ppi.neighbors(v)):
            if node_to_cluster_map.__contains__(n) and v[0:2] == n[0:2]:
                if not related_clusters.__contains__(node_to_cluster_map[n]) and not cluster.__contains__(n):
                    related_clusters.append(node_to_cluster_map[n])
    internal_conserved /= 2
    for i in related_clusters:
        result = add_tuples(result, ciq_between_two(cluster, clusters[i]))
    return result, internal_conserved

def add_tuples(t1, t2):
    result = (t1[0] + t2[0], t1[1] + t2[1])
    return result

def ciq_between_two(c1, c2):
    mutual_k = 0
    n = 0
    kset_t = []
    for nn in c1:
        kname = nn[0:2]
        if not kset_t.__contains__(kname):
            kset_t.append(kname)
    kset2 = []
    for nn in c2:
        if not kset2.__contains__(nn[0:2]) and kset_t.__contains__(nn[0:2]):
            kset2.append(kname)
    mutual_k = len(kset2)
    kset_t.clear()
    # calculate edge weights
    for nn1 in c1:
        for nn2 in c2:
            if ppi.has_edge(nn1, nn2) and nn1[0:2] == nn2[0:2]:
                if not kset_t.__contains__(nn1[0:2]):
                    kset_t.append(nn1[0:2])
                n = n + 1
    e = len(kset_t)
    weight = 0.0
    if e < 2:
        weight = 0.0
    elif mutual_k != 0:
        weight = (e * 1.0) / (mutual_k * 1.0)
    else:
        weight = 0.0
    result = (n * weight, n)
    return result

if __name__ == '__main__':

    starttime = time.time()

    datafile = "data-real.txt"
    global beta
    global alpha
    global ppi
    global sim_ppi
    global node_prio
    global simmap
    global node_to_cluster_map
    global NODEs_clique
    global NODEs_Weight
    global candidate_clique

    global matchset
    global aligned_node
    global labelnode
    global cluster_node_map
    global clusters
    global scores
    for i in range(0,11):
        alpha = float(i / 10)
        print(alpha)
        for i in range(0,11):
            ppi = nx.Graph()
            sim_ppi = nx.Graph()
            node_prio = []
            simmap = {}
            node_to_cluster_map = {}
            NODEs_clique = {}
            NODEs_Weight = {}
            candidate_clique = []

            matchset = []
            aligned_node = []
            labelnode = []
            cluster_node_map = {}
            clusters = {}
            scores = {}
            beta = float(i / 10)
            print(beta)
            resultfile = "res-" + str(alpha)+"-"+str(beta)+".txt"
            print("alpha is {alpha}, beta is {beta}".format(alpha=alpha,beta=beta))
            readnetwork(datafile)
            # PPI network and score reading finished
            print("start clique...")
            createclique()
            print("find clique done!")
            # clique_weight = sorted(NODEs_Weight.items(), key = lambda kv:(kv[1], kv[0]))

            # for i in range(0,int(len(clique_weight)/2)):
            #     key = clique_weight.pop(0)
            #     print(key)
            #     clique = NODEs_clique.pop(key[0])
            #     if clique:
            #         candidate_clique.remove(clique)
            simuldated_annealing(candidate_clique)

            fw = open(resultfile, 'w')
            for cluster in clusters.values():
                cluster.sort()
                for v in cluster:
                    fw.write(v + ' ')
                fw.write('\n')
            endtime = time.time()
            dtime = endtime - starttime
            print(dtime)

