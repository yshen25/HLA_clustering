#! usr/bin/env python3
"""
assess cluster-wise stability of hierarchical clustering results
using boot striping and Jaccard index
Using linear balanced assignment problem algorithm, a variation of stable marriage problem
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

def jaccard(member_index1:list, member_index2:list) -> float:
    # calculate jaccard index between two groups
    s1 = set(member_index1)
    s2 = set(member_index2)
    return len(s1.intersection(s2)) / len(s1.union(s2))

def min_jaccard(clustering1:np.array, clustering2:np.array):
    groups1 = np.unique(clustering1) # name of groups
    groups2 = np.unique(clustering2)
    jar_matrix = np.zeros((len(groups1),len(groups2))) # jaccard index matrix, storing jaccard index of all-to-all groups

    for i in groups1:
        for j in groups2:
            members1 = np.where(clustering1 == i) # index of group members
            members2 = np.where(clustering2 == j)
            jar_matrix[i,j] = jaccard(members1, members2)

    row_ind, col_ind = linear_sum_assignment(jar_matrix) # find the solution with minimum sum
    return jar_matrix[row_ind, col_ind].sum()