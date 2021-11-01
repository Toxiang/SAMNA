import networkx as nx
import numpy as np
import random
import math

class SA:
    current_score = 0.0

sim_ppi = nx.Graph()
matchset = []
aligned_node = []
labelnode = []
cluster_node_map = {}
clusters = []
simmap = {}

def simuldated_annealing(candidate):
    t_Max = 100  #initiate temperature
    t_Min = 10   #minimum value of terperature
    k_Max = 100
    k = 100   #number of iterations
    s = 0.005
    n_Max = 100

    k = 0
    t = t_Max
    step = (t_Max - t_Min) / k_Max

    while k <= k_Max:
        t = t - step
        beta = 1.0 / (s * t)
        for i in n_Max:
            delta = 0.0
            select = random.choice(candidate)
            overlapValue, delta, clu = isOverlapped(select)
            if not overlapValue:
                add(select,clu)
            else:
                if delta > 0 or random.random()<math.exp(beta*delta):
                    updateState(select,delta,overlapValue)


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
                    newmp2.move(node)
            deltaScore[1] = getDeltaScore(newmp1, cluster1, c_n1) + getDeltaScore(newmp2, cluster2,c_n2)
            clu_index[1] = c_n1
            index_i = 1
        else:
            for u in newmp1:
                if u[0:2]==species2:
                    newmp1 = newmp1.replace(u, protein2)
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
                    newmp2 = newmp2.replace(u, protein1)
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
                    newmp1 = newmp1.replace(u, protein2)
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
                    newmp2 = newmp2.replace(u, protein1)
                    deltaScore[8] = getDeltaScore(newmp2, cluster2,c_n2)
                    clu_index[8] = c_n2
                    overlapValue = 8
    myDeltaScore = deltaScore[overlapValue]
    cluster_index = clu_index[overlapValue]
    return overlapValue, myDeltaScore,cluster_index

def isCombine(matchset, protein, species):
    for u in matchset:
        if u[0:2] == species:
            if not sim_ppi.has_edge(u,protein):
                return False
        else:
            return True
    return True


def add(candidate,i):
    protein1 = candidate[0]
    protein2 = candidate[1]
    labelnode.append(protein1)
    labelnode.append(protein2)
    cluster_node_map[protein1] = i
    cluster_node_map[protein2] = i
    if protein1 not in clusters[i]:
        clusters[i].append(protein1)
    if protein2 not in clusters[i]:
        clusters[i].append(protein2)

def getDeltaScore(newcluster,oldcluster,i):
    score1 = 0.0
    score2 = 0.0
    icq_new, kk = calculate_icq(newcluster)
    int_cons = 0
    ciq_new, int_cons = calculate_ciq(newcluster)



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


def updateState(candidate,delta,overlapValue):
    protein1 = candidate[0]
    protein2 = candidate[1]
    if overlapValue == 1:
        species2 = protein2[0:2]
