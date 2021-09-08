#!usr/bin/env python3
# -*- coding: utf-8 -*-

from sklearn.cluster import SpectralClustering, DBSCAN, AffinityPropagation
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

def cluster(InCSV):
    Mat = pd.read_csv(InCSV, index_col=0)
    Mat = Mat.add(Mat.T, fill_value=0)
    clustering = DBSCAN(eps=9, min_samples=2, metric="precomputed").fit(Mat)
    print(clustering.labels_)

    return

def dendro(InCSV, OutTree):
    Mat = pd.read_csv(InCSV, index_col=0)
    Mat = Mat.add(Mat.T, fill_value=0)
    dists = squareform(Mat)
    # print(dists)
    A = linkage(dists, "single")
    # print(A)
    tree = to_tree(A,False)
    OutFile = getNewick(tree, "", tree.dist, Mat.index)
    with open(OutTree, "w") as fh:
        fh.write(OutFile)
    # dendrogram(A, labels=Mat.index)
    # plt.title("HLAA")
    # plt.rcParams['axes.labelsize'] = 4
    # plt.show()
    return

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

# cluster("HLAA_distance.csv")
# dendro("HLAA_distance.csv", "HLAA_tree.newick")
dendro("HLAA_distance_s.csv", "HLAA_tree_s.newick")
dendro("HLAB_distance_s.csv", "HLAB_tree_s.newick")