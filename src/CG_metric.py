#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate distances between coarse-grained structures
"""

import os
import glob
from pathlib import Path
from itertools import combinations, combinations_with_replacement
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from .PropertyParams import GranthamSim, SM_THREAD_NORM_Sim, PMBEC_Sim

class CGCalculator():

    def __init__(self, SimilarityMatrix="Grantham") -> None:
        
        # # fixed parameters
        self.sigma = 0.5
        self.k = 1

        if SimilarityMatrix == "Grantham":
            self.SimMtx = GranthamSim

        elif SimilarityMatrix == "SM_THREAD_NORM":
            self.SimMtx = SM_THREAD_NORM_Sim

        elif SimilarityMatrix == "PMBEC":
            self.SimMtx = PMBEC_Sim

        else:
            raise ValueError(f"Similarity Matrix not recognized: {SimilarityMatrix}")
            
    def ParamExtract(self, DATFile):
        DAT = pd.read_csv(DATFile)

        resi = DAT['Residue'].values.reshape((-1,1))
        
        # convert back HIS variants
        # resi = resi.reshape(-1).tolist()
        resi = ["HIS" if s in ["HID", "HIE", "HIP"] else s for s in resi]
        resi = np.array(resi)

        coord = DAT[['X', 'Y', 'Z']].values
        weight = DAT['Weight'].values.reshape((-1,1))

        accepted = np.flatnonzero(weight) # remove unaccepted residues with weight 0

        return coord[accepted], resi[accepted].reshape(-1).tolist(), weight[accepted]

    def ResiPairSim(self, ResiA, ResiB):
        """
        Residue similarity
        """
        ResiPairComb = [(x, y) for x in ResiA for y in ResiB]
        # print(ResiPairComb)
        ResiPairSim = np.array([self.SimMtx[i][j] for i,j in ResiPairComb])

        return ResiPairSim.reshape((len(ResiA), len(ResiB)))

    def CloudSimilarity(self, CoordA, ResiA, CoordB, ResiB, WeightA, WeightB):
        
        ResiPairSim_score = self.ResiPairSim(ResiA, ResiB)

        SimScore = np.sum( np.reciprocal(np.cosh(self.sigma*cdist(CoordA, CoordB, "euclidean")))**self.k * ResiPairSim_score * np.outer(WeightA, WeightB) )

        return SimScore

    def CalcSim(self, DATpair):
        
        """Similarity score between two point clouds"""

        # print(f"Simi: {comb}")
        
        CoordA, ResiA, WeightA = self.ParamExtract(DATpair[0])
        CoordB, ResiB, WeightB = self.ParamExtract(DATpair[1])
        
        return (DATpair, self.CloudSimilarity(CoordA, ResiA, CoordB, ResiB, WeightA, WeightB))
    
    def SaveDist(self, OutCSV):
        
        self.DistMat.to_csv(OutCSV)
        
        return
    
    def CalcDist(self, DATDir, ListFile=None):
        """
        Distance of two molecules in PDB files
        ======================================
        Input:
            DATDir: Directory of coarse-grained HLA structures
            ListFile (optional): List file for selecting alleles, see "../Dataset_split" directory
        """
        DATList = [InDAT for InDAT in glob.glob(f"{DATDir}/*.csv")]
        if ListFile:
            with open(ListFile, "r") as ALF:
                RequiredList = [f"{DATDir}/{line.strip()}" for line in ALF]
                if set(RequiredList)-set(DATList):
                    raise ValueError(f"File in list is not found: {list(set(RequiredList)-set(DATList))}")
                else:
                    DATList = list(set(RequiredList)&set(DATList))
        lite_DATList = [Path(DAT).stem for DAT in DATList] # remove directory name and suffix
        AlleleComb_wi = combinations_with_replacement(DATList, 2)
        AlleleComb_wo = combinations(DATList, 2)

        self.DistMat = pd.DataFrame(np.zeros((len(lite_DATList), len(lite_DATList))), index=lite_DATList, columns=lite_DATList)

        SimilarityMat = {}

        pool = Pool(os.cpu_count())
        result = pool.map(self.CalcSim, AlleleComb_wi)
        
        pool.close()
        pool.join()
        
        for i in result:
            SimilarityMat[i[0]] = i[1]
        
        for comb in AlleleComb_wo:
            # print(f"Dist: {comb}")
            distance = np.sqrt(SimilarityMat[(comb[0], comb[0])] + SimilarityMat[(comb[1], comb[1])] - 2 * SimilarityMat[comb])
            self.DistMat.loc[Path(comb[1]).stem, Path(comb[0]).stem] = distance

        return