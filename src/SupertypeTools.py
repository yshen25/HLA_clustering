#!usr/bin/env python3
"""
High-level tools for HLA clustering
"""
from sklearn.cluster import DBSCAN, AgglomerativeClustering
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, to_tree

import matplotlib.pyplot as plt
import seaborn as sn
from scipy.stats import linregress

from .CalcScore import Calculator
from .CG_metric import CGCalculator

from Bio import SeqIO, Align
from itertools import combinations

# from sklearn.linear_model import LinearRegression

# ==== Calculating distance matrix ====
def CalcMat(DATDir, AlleleListFile, contact, weight):
    """
    Using all-atom distance metric
    """

    calc = Calculator(DATDir, AlleleListFile, ContactResi=contact, ResiWeight=weight)

    calc.CalcDist()

    return calc.DistMat

def CGCalcMat(CGDATDir, AlleleListFile, SimMtx="Grantham", sigma=None, k=None, DistMat_output=None, Standardize:bool=False):
    """
    Using coarse-grained distance metric
    ====================================
    Input:
        CGDATDir: Directory of coarse-grained HLA structures
        SimMtx (optional): Similarity matrix, choose from ["Grantham", "SM_THREAD_NORM", "PMBEC"]
        AlleleListFile (optional): List file for selecting alleles, see "../Dataset_split" directory
        sigma, k (optional): Shape parameters
        DistMat_output (optional): File name of distance matrix
        Standardize (optional): If true, standardize the distance matrix to [0-1]

    Output:
        Distance_matrix
    """
    metric = CGCalculator(SimilarityMatrix=SimMtx)
    if sigma:
        metric.sigma = sigma

    if k:
        metric.k = k

    metric.CalcDist(CGDATDir, AlleleListFile)

    if Standardize:
        metric.DistMat = (metric.DistMat - metric.DistMat.min().min()) / (metric.DistMat.max().max() - metric.DistMat.min().min())

    if DistMat_output:
        metric.SaveDist(DistMat_output)

    return metric.DistMat

def MSAMat(InFile, scale=1):
    """
    Similarity from pairwise sequence alignment
    """
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

# ==== Clustering ====
def hierarchical_cluster(Mat, N=None, square=False, L='complete', threshold=None, outtree=None, plot_dendro=False, centroid=False, color_threshold=None):
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    model = AgglomerativeClustering(n_clusters=N, affinity='precomputed', linkage=L, distance_threshold=threshold, compute_distances=True, compute_full_tree=True).fit(Mat)
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

def DBSCAN_cluster(Mat:pd.DataFrame, epsilon:float, square=False, MinSample=5):
    # Mat = pd.read_csv(InCSV, index_col=0)
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    clustering = DBSCAN(eps=epsilon, min_samples=MinSample, metric="precomputed", n_jobs=-1).fit(Mat)
    labels = clustering.labels_
    result = pd.Series(labels, index=Mat.index)

    return result

# ==== Visualization tools ====
def dist_heatmap(Mat, square=True, order=None, size=(10,10), label=False, line=False, labelsize=8, **cbar_kw):
    """
    Visualize distance matrix as heatmap
    """

    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)

    if order:
        if type(order[0]) == list:
            flat_order = [item for sublist in order for item in sublist]
        else:
            flat_order = order
        Mat = Mat[flat_order] # re-arrange row order
        Mat = Mat.reindex(flat_order) # re-arrange column order

    # standardize to [0,1]
    Mat = (Mat - Mat.min().min()) / (Mat.max().max() - Mat.min().min())
    
    plt.figure(figsize=size)
    if label:
        ticks = True
    else:
        ticks = False

    g = sn.heatmap(Mat, square=True, xticklabels=ticks, yticklabels=ticks, cbar_kws=cbar_kw)
    g.axes.tick_params(axis='both', labelsize=labelsize, pad=50)

    # seperate lines between
    if line:
        split = np.cumsum([len(sublist) for sublist in order])
        for line in split[:-1]:
            plt.axhline(y=line, color='k', linestyle='-')
            plt.axvline(x=line, color='k', linestyle='-')
    plt.show()
    return


def getNewick(node, newick, parentdist, leaf_names):
    """
    convert dendrogram to newick format
    """
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
    """
    Create linkage matrix and then plot the dendrogram
    """
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

    # Plot the corresponding dendrogram
    plt.figure(figsize=(80,10))

    dendro = dendrogram(linkage_matrix, truncate_mode=truncate, labels=labels, leaf_font_size=12, get_leaves=True, color_threshold=color_threshold)

    plt.show()

    if outtree:
        tree = to_tree(linkage_matrix)
        OutFile = getNewick(tree, "", tree.dist, labels)
        with open(outtree, "w") as fh:
            fh.write(OutFile)

    return dendro['leaves']

def square2triangle(InMat):
    keep = np.invert(np.triu(np.ones(InMat.shape)).astype('bool'))
    return InMat.mask(keep, other=0)

def triangle2square(InMat):
    return InMat.add(InMat.T, fill_value=0)

def crop_mtx(Mtx:pd.DataFrame, order:list, flatten:bool=False):
    
    if type(order[0]) == list:
        flat_order = [item for sublist in order for item in sublist]
    else:
        flat_order = order
    
    Mtx = Mtx.loc[flat_order, flat_order] # re-arrange row order

    if flatten:
        # since matrix is symmetric, each pairwise distance exist 2 copies
        # first extract lower triangle (no diagonal), then extract values to make sure only extract unique values
        keep = np.invert(np.triu(np.ones(Mtx.shape)).astype('bool')).flatten()
        return Mtx.to_numpy().flatten()[keep]

    return Mtx

def correlation(ArrayA, ArrayB, show_plot=True):
    """
    Correlation plot between two array
    =============================
    Input:
        ArrayA, ArrayB: 1-D array with same shape
        show_plot: if true, draw correlation plot
    return:
        (slope, intercept, rvalue)
    """

    # standardize to 0-1
    xx = (ArrayA - np.min(ArrayA)) / (np.max(ArrayA) - np.min(ArrayA))
    yy = (ArrayB - np.min(ArrayB)) / (np.max(ArrayB) - np.min(ArrayB))

    # xx = xx.reshape(-1, 1)
    # yy = yy.reshape(-1, 1)

    # rgs = LinearRegression(fit_intercept=True)
    # rgs.fit(xx, yy)
    # slope, intercept, rvalue = rgs.coef_[0][0], rgs.intercept_[0], rgs.score(xx, yy)    

    slope, intercept, rvalue, _, _ = linregress(xx, yy)

    # equation = f"Y = {round(slope, 3)}X"

    if intercept >= 0:
        equation = f"Y = {round(slope, 3)}X+{round(intercept, 3)}"
    else:
        equation = f"Y = {round(slope, 3)}X{round(intercept, 3)}"
    # intercept = 0
    if show_plot:
        plt.figure(figsize=(10,10))
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.scatter(xx, yy)

        x_vals = np.array([0, 1])
        y_vals = intercept + slope * x_vals
        # y_vals = slope * x_vals
        plt.plot(x_vals, y_vals, '--')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel("Structure distance", fontsize=24)
        plt.ylabel("Peptide binding specificity distance", fontsize=24)

        plt.text(0.02,0.95, equation, fontsize=24)
        plt.text(0.02,0.9, f"R = {round(rvalue, 3)}", fontsize=24)
        plt.show()
    
    return (slope, intercept, rvalue)

def SSE(Mat:pd.DataFrame, clusters:list, square=True):
    # SSE divided by number of elements in each cluster
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    
    sum_SSE = 0

    for clst in clusters:
        
        clst_mat = Mat.loc[clst, clst]
        centroid = clst_mat.sum().idxmin()
        sum_SSE += clst_mat[centroid].pow(2).sum()
    
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

def Tuning_N(ClusterMat, Nmin, Nmax, RefMat=None, ClusterSilhouette=False, RefSilhouette=False) -> tuple:
    """
    Determining optimized number of clusters
    Arguments:
        ClusterMat: Distance matrix that the clustering is based on, lower triangular form
        RefMat: Distance matrix that the performance accessment is based on, lower triangular form
        Nmin: starting number of clusters
        Nmax: ending number of clusters
    Return:
        (StructSSE, BASSE, StructSilhouette, BASilhouette)
    """
    ClusterSSE = []
    RefSSE = []
    Silhouette_C = []
    Silhouette_R = []

    # dist_list = []

    for i in range(Nmin, Nmax+1):

        # initialize optional parameters
        # BA_err = 'NA'
        # SilhouetteScore = 'NA'
        # BASilhouetteScore = 'NA'

        cluster, _ = hierarchical_cluster(ClusterMat, square=True, N=i, L='complete', threshold=None)
        #complete average single
        groups = [group[1].index.tolist() for group in cluster.groupby(cluster)]
        # print(groups)
        
        ClusterSSE.append(SSE(ClusterMat, groups))

        if RefMat is not None:
            RefSSE.append(SSE(RefMat, groups))

        if ClusterSilhouette:
            Silhouette_C.append(Silhouette(ClusterMat, groups))

        if RefSilhouette:
            Silhouette_R.append(Silhouette(RefMat, groups))

    return ClusterSSE, RefSSE, Silhouette_C, Silhouette_R

def elbow_plot(Struct_Mat, BA_Mat, Additional_Bar_group:list=None, Nmin=1, Nmax=12):
    """
    Draw elbow plot (Sum-of-Squared-Error (SSE) versus number of clusters (N)) to determine number of clusters
    """
    
    SSE_struct, SSE_BA = Tuning_N(Struct_Mat, BA_Mat, Nmin, Nmax)
    xx = range(Nmin, Nmax+1)

    lines = []
    labels = []

    fig, ax1 = plt.subplots(figsize=(6,10))
    ax2 = ax1.twinx()

    line1, = ax1.plot(xx, SSE_BA, c='b', marker='^', mfc='None', mec='b', ms='8', mew=3, alpha=0.6, label="NetMHCpan SSE")
    line2, = ax2.plot(xx, SSE_struct, c='g', marker='v', mfc='None', mec='g', ms='8', mew=3, alpha=0.6, label="Structure SSE")
    lines.append(line1)
    labels.append("NetMHCpan SSE")
    lines.append(line2)
    labels.append("Structure SSE")

    ax1.set_xlabel('Number of clusters', fontsize=20)
    ax1.set_xticks(range(1,Nmax+1,2))
    ax1.tick_params(axis='x', labelsize=16)
    
    ax1.set_ylabel('Binding peptide SSE', color='tab:blue', fontsize=20)
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=16)

    ax2.set_ylabel('Structure distance SSE', color='tab:green', fontsize=20)
    ax2.tick_params(axis='y', labelcolor='tab:green', labelsize=16)

    if Additional_Bar_group:
        for group in Additional_Bar_group:
            lines.append(ax1.bar(group[0], group[1], alpha=0.6, label=group[2], linewidth=2, edgecolor='b'))
            labels.append(f"{group[2]} (n={group[0]})")

    ax1.legend(lines, labels, prop={"size":16})
    ax1.grid(linestyle='--')

    # fig.legend()
    plt.show()

    return