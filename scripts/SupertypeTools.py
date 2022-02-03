from turtle import right
from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sn

from skbio import DistanceMatrix
from skbio.tree import nj

from CalcScore import Calculator
from CGmodel import CGCalculator

from Bio import SeqIO, Align
from itertools import combinations
# ==== full atom tools ====
def CalcMat(DATDir, AlleleListFile, OutCSV, contact, weight):
    # l: charge param
    # s: spatial param
    # d: depth param
    calc = Calculator(DATDir, AlleleListFile, OutCSV, ContactResi=contact, ResiWeight=weight)
    # calc.l = l
    # calc.sigma = s
    # calc.d = d
    calc.CalcDist()

    return calc.DistMat

# ==== coarse grain tools ====
def CGCalcMat(DATDir, AlleleListFile, OutCSV, contact, weight, pairwise=False):

    calc = CGCalculator(DATDir ,AlleleListFile, OutCSV, ContactResi=contact, ResiWeight=weight, Pairwise=pairwise)

    calc.CalcDist()

    return calc.DistMat

# ==== universal tools ====
def heatmap(Mat, order=None, size=(10,10), label=False, line=False):
    # Mat = pd.read_csv(InCSV, index_col=0)
    Mat = Mat.add(Mat.T, fill_value=0)
    # print(Mat.index)

    if order:
        if line:
            flat_order = [item for sublist in order for item in sublist]
        else:
            flat_order = order
        Mat = Mat[flat_order] # re-arrange row order
        Mat = Mat.reindex(flat_order) # re-arrange column order
    # print(Mat.index, Mat.columns)

    # fig, axs = plt.subplots(figsize=size)
    plt.figure(figsize=size)
    
    g = sn.heatmap(Mat, square=True, xticklabels=True, yticklabels=True, cbar_kws={"shrink": .8})
    g.axes.tick_params(axis='both', labelsize=8, pad=45)
    if label:
        g.axes.set_xticklabels(labels=label,va='bottom',ha='center')
        g.axes.set_yticklabels(labels=label,va='center',ha='left')

    else:
        g.axes.set_xticklabels(labels=g.axes.get_xticklabels(),va='bottom',ha='center')
        g.axes.set_yticklabels(labels=g.axes.get_yticklabels(),va='center',ha='left')

    # seperate lines between
    if line:
        split = np.cumsum([len(sublist) for sublist in order])
        # print(split)
        for line in split[:-1]:
            plt.axhline(y=line, color='k', linestyle='-')
            plt.axvline(x=line, color='k', linestyle='-')
    plt.show()
    return

def DBSCAN_cluster(InCSV):
    Mat = pd.read_csv(InCSV, index_col=0)
    Mat = Mat.add(Mat.T, fill_value=0)
    clustering = DBSCAN(eps=9, min_samples=2, metric="precomputed").fit(Mat)
    print(clustering.labels_)

    return

"""
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

def dendro(Mat, OutTreeFile):
    # Mat = pd.read_csv(InCSV, index_col=0)
    Mat = Mat.add(Mat.T, fill_value=0)
    dists = squareform(Mat)
    # print(dists)
    A = linkage(dists, "single")
    # print(A)
    tree = to_tree(A,False)
    OutFile = getNewick(tree, "", tree.dist, Mat.index)
    with open(OutTreeFile, "w") as fh:
        fh.write(OutFile)
    # dendrogram(A, labels=Mat.index)
    # plt.title("HLAA")
    # plt.rcParams['axes.labelsize'] = 4
    # plt.show()
    return
"""

def Matrix2Dendro(Mat, OutTreeFile=None, label=None):
    Mat = Mat.add(Mat.T, fill_value=0)
    
    # if label:
    #     Mat = pd.DataFrame(Mat.values, index=label, columns=label)
    #     print(Mat)

    dm = DistanceMatrix(Mat, Mat.columns)
    tree = nj(dm)
    # print(tree.ascii_art())

    if OutTreeFile:
        result = str(tree)

        if label:
            for old_string, dst_string in zip(Mat.columns, label):
                result = result.replace(old_string, dst_string)
        
        with open(OutTreeFile, "w") as out:
            out.write(result)
    
    return

def MSAMat(InFile, scale=0.01):
    allele_dict = SeqIO.to_dict(SeqIO.parse(InFile, "fasta"))

    # AlleleComb_wi = combinations_with_replacement(DATList, 2)
    allele_list = list(allele_dict.keys())
    AlleleComb_wo = combinations(allele_list, 2)
    DistMat = pd.DataFrame(np.zeros((len(allele_list), len(allele_list))), index=allele_list, columns=allele_list)
    # print(AlleleComb_wo)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    for comb in AlleleComb_wo:
        DistMat.loc[comb[1],comb[0]] = np.reciprocal(aligner.score(allele_dict[comb[1]].seq, allele_dict[comb[0]].seq)*scale)

    return DistMat