# from turtle import right
from sklearn.cluster import DBSCAN, AgglomerativeClustering
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, to_tree
# from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sn

from skbio import DistanceMatrix
from skbio.tree import nj

from CalcScore import Calculator
from CGmodel import CGCalculator

from Bio import SeqIO, Align
from itertools import combinations, count
# ==== full atom tools ====
def CalcMat(DATDir, AlleleListFile, contact, weight):
    # l: charge param
    # s: spatial param
    # d: depth param
    calc = Calculator(DATDir, AlleleListFile, ContactResi=contact, ResiWeight=weight)
    # calc.l = l
    # calc.sigma = s
    # calc.d = d
    calc.CalcDist()

    return calc.DistMat

# ==== coarse grain tools ====
def CGCalcMat(DATDir, AlleleListFile, contact, weight, sigma=None, w=None, pairwise=False):

    calc = CGCalculator(DATDir ,AlleleListFile, ContactResi=contact, ResiWeight=weight, Pairwise=pairwise)
    if sigma:
        calc.sigma = sigma

    if w:
        calc.w = w

    calc.CalcDist()

    return calc.DistMat

# ==== universal tools ====
def heatmap(Mat, square=False, order=None, size=(10,10), label=False, line=False, labelsize=8, **cbar_kw):
    sn.set(font_scale=2)
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    # print(Mat.index)

    if order:
        if type(order[0]) == list:
            flat_order = [item for sublist in order for item in sublist]
        else:
            flat_order = order
        Mat = Mat[flat_order] # re-arrange row order
        Mat = Mat.reindex(flat_order) # re-arrange column order
    # print(Mat.index, Mat.columns)

    # fig, axs = plt.subplots(figsize=size)
    plt.figure(figsize=size)
    if label:
        ticks = True
    else:
        ticks = False

    g = sn.heatmap(Mat, square=True, xticklabels=ticks, yticklabels=ticks, cbar_kws=cbar_kw)
    g.axes.tick_params(axis='both', labelsize=labelsize, pad=50)
    if label:
        g.axes.set_xticklabels(labels=label,va='bottom',ha='center')
        g.axes.set_yticklabels(labels=label,va='center',ha='left')

    # else:
    #     g.axes.set_xticklabels(labels=g.axes.get_xticklabels(),va='bottom',ha='center')
    #     g.axes.set_yticklabels(labels=g.axes.get_yticklabels(),va='center',ha='left')

    # seperate lines between
    if line:
        split = np.cumsum([len(sublist) for sublist in order])
        # print(split)
        for line in split[:-1]:
            plt.axhline(y=line, color='k', linestyle='-')
            plt.axvline(x=line, color='k', linestyle='-')
    plt.show()
    return

def DBSCAN_cluster(Mat:pd.DataFrame, epsilon:float, square=False, MinSample=5):
    # Mat = pd.read_csv(InCSV, index_col=0)
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    clustering = DBSCAN(eps=epsilon, min_samples=MinSample, metric="precomputed", n_jobs=-1).fit(Mat)
    labels = clustering.labels_
    result = pd.Series(labels, index=Mat.index)

    return result

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

def plot_dendrogram(model, truncate, labels, color_threshold, outtree=None):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    # linkage_matrix = linkage(Mat, method='average')

    # Plot the corresponding dendrogram
    plt.figure(figsize=(80,10))
    # fig, ax = plt.subplots(figsize=(80,10))
    dendro = dendrogram(linkage_matrix, truncate_mode=truncate, labels=labels, leaf_font_size=12, get_leaves=True, color_threshold=color_threshold)

    plt.show()

    if outtree:
        tree = to_tree(linkage_matrix)
        OutFile = getNewick(tree, "", tree.dist, labels)
        with open(outtree, "w") as fh:
            fh.write(OutFile)

    return dendro['leaves']

def hierarchical_cluster(Mat, N=None, square=False, L='complete', threshold=None, outtree=None, plot_dendro=False, centroid=False, color_threshold=None):
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    model = AgglomerativeClustering(n_clusters=N, affinity='precomputed', linkage=L, distance_threshold=threshold, compute_distances=True).fit(Mat)
    result = pd.Series(model.labels_, index=Mat.index)

    order = None
    
    if plot_dendro:
        order = plot_dendrogram(model, None, Mat.index, color_threshold, outtree)

    if centroid:
        centers = []
        for i in result.groupby(by=result):
            group = i[1].index.to_numpy()
            centers.append(Mat.loc[group,group].sum(axis=0).idxmin())

        return result, order, pd.Series(range(len(centers)), index=centers)

    return result, order

def robust_cluster():
    return


"""
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

def Matrix2Dendro(Mat, square=False, OutTreeFile=None, label=None):
        
    if not square:
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

def MSAMat(InFile, scale=1):
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

def SSE(mat:pd.DataFrame, groups:list):
    # mat = pd.read_csv(MatrixCSV, index_col=0)
    # SSE divided by number of elements in each cluster
    sum_SSE = 0

    for group in groups:
        # print(group)
        group_square = mat.loc[group,group]
        sum_SSE += group_square.values.sum()/group_square.shape[0]
    
    return sum_SSE

def Silhouette(Mat:pd.DataFrame, groups:list, square=False):
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    score = []

    if len(groups) == 1:
        return 0

    for i in range(len(groups)):

        if len(groups[i]) == 1:
            # if only one element in a group, the silhouette score is 0 (arbitrary)
            score.append(0)
            continue

        out_groups = groups[0:i] + groups[i+1:]
        # out_groups.remove(group)

        for allele in groups[i]:
            
            # distance within groups
            ai = Mat.loc[groups[i],allele].sum() / (len(groups[i]) - 1)
            # distance to neighbor cluster
            bi = np.min([Mat.loc[out_group,allele].sum() / len(out_group) for out_group in out_groups])

            if ai <= bi:
                score.append(1-ai/bi)

            else:
                score.append(bi/ai-1)

        # print(score)

    return np.mean(score)