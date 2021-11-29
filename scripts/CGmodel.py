#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From full-atom DAT file to coarse-grained DAT file
Residues are represented by side-chain center-of-mass
"""
import os
import time
from scipy.spatial.distance import cdist
from itertools import combinations, combinations_with_replacement
from multiprocessing import Pool

import numpy as np
# from scipy.stats import hypsecant
import pandas as pd

from PropertyParams import AtomicMass, AASim

def CenterOfMass(DATFile, OutFile):
    """
    center-of-mass of side chain
    """
    # if OutFile is None:
    #     OutFile = DATFile.split(".")[0] + "_CG.csv"

    DAT = pd.read_csv(DATFile)
    DAT = DAT[~DAT.Atom.isin(["N", "CA", "C", "O"])] # filter out backbone atoms
    ResGen = DAT.groupby(by='ResNum')

    CG_list = []
    for resi in ResGen:
        x = y = z = 0
        total_mass = 0
        resnum = resi[0]
        resatoms = resi[1]
        # print(resatoms)
        resname = resatoms["Residue"].iloc[0]

        for atom in resatoms[["Atom", "X", "Y", "Z"]].to_numpy():
            mass = AtomicMass[atom[0][0]] # real atom name is the first letter of pdb atom name
            total_mass += mass
            x += mass * atom[1]
            y += mass * atom[2]
            z += mass * atom[3]


        CG_list.append([resnum, resname, x/total_mass, y/total_mass, z/total_mass])

    CG_DAT = pd.DataFrame(CG_list, columns=["ResNum", "Residue", "X", "Y", "Z"])
    CG_DAT.to_csv(OutFile, index=False)

    return

# CenterOfMass("crystal/A_mean/DAT/A01_01.csv")

class CGCalculator():
    """
    Calculate distance between CG models
    """
    def __init__(self, DATDir, OutCSV, ContactResi:list=None, ResiWeight:dict=None) -> None:
        self.OutCSV = OutCSV
        self.DATDir = DATDir
        self.AASimDict = AASim

        self.ContactResi = ContactResi
        self.Weight = ResiWeight
        self.OnlyGroove = 0

        # Allele combinations
        DATList = [a.split(".")[0] for a in sorted(os.listdir(DATDir))]
        self.AlleleComb_wi = combinations_with_replacement(DATList, 2)
        self.AlleleComb_wo = combinations(DATList, 2)

        # Distance matrix for output, in lower triangular form
        self.DistMat = pd.DataFrame(np.zeros((len(DATList), len(DATList))), index=DATList, columns=DATList)

    def ResiPairSim(self, ResiA, ResiB):
        ResiPairComb = [(x, y) for x in ResiA for y in ResiB]
        ResiPairSim = np.array([AASim[i][j] for i,j in ResiPairComb])

        return ResiPairSim.reshape((len(ResiA), len(ResiB)))


    def CloudSimilarity(self, CoordA, ResiA, CoordB, ResiB, WeightA, WeightB):
        
        ResiPairSim_score = self.ResiPairSim(ResiA, ResiB)

        SimScore = np.sum( np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))) * ResiPairSim_score * np.outer(WeightA, WeightB) )

        return SimScore

    def ParamExtract(self, DATFile):
        DAT = pd.read_csv(DATFile)
        resnum = DAT['ResNum'].values.reshape((-1,1))
        resi = DAT['Residue'].values.reshape((-1,1))
        coord = DAT[['X', 'Y', 'Z']].values

        return resnum, resi, coord

    def AssignWeight(self, resnum, WeightDict):
        weight = []
        for i in resnum.reshape(-1):
            if i in WeightDict:
                weight.append(WeightDict[i])

            else:
                weight.append(1)

        return weight

    def CalcSim(self, comb):
        
        """Similarity score between two point clouds"""

        # print(f"Simi: {comb}")
        
        resnumA, resiA, CoordA = self.ParamExtract(f"{self.DATDir}/{comb[0]}.csv")
        resnumB, resiB, CoordB = self.ParamExtract(f"{self.DATDir}/{comb[1]}.csv")
        # print(CoordA.shape)

        # filter atoms
        if self.ContactResi:
            passA = np.isin(resnumA, self.ContactResi).reshape(-1)
            resnumA = resnumA[passA]
            resiA = resiA[passA]
            CoordA = CoordA[passA]

            passB = np.isin(resnumB, self.ContactResi).reshape(-1)
            resnumB = resnumB[passB]
            resiB = resiB[passB]
            CoordB = CoordB[passB]

        # weigh atoms
        if self.Weight:

            # invert keys and values
            WeightDict = {}
            for key, values in self.Weight.items():
                for v in values:
                    WeightDict[v] = key
            WeightA = self.AssignWeight(resnumA, WeightDict)
            WeightB = self.AssignWeight(resnumB, WeightDict)

        else:
            WeightA = np.ones_like(resnumA)
            WeightB = np.ones_like(resnumB)

        # convert numpy ndarray of objects to plain list
        resiA = resiA.reshape(-1).tolist()
        resiB = resiB.reshape(-1).tolist()

        # convert back HIS variants
        resiA = ["HIS" if s in ["HID", "HIE", "HIP"] else s for s in resiA]
        resiB = ["HIS" if s in ["HID", "HIE", "HIP"] else s for s in resiB]
        
        return (comb, self.CloudSimilarity(CoordA, resiA, CoordB, resiB, WeightA, WeightB))

    def CalcDist(self):
        
        """
        Distance of two molecules in PDB files
        """
        
        # Similarity score of two alleles
        SimilarityMat = {}

        pool = Pool(os.cpu_count())
        # pool = Pool(1)
        result = pool.map(self.CalcSim, self.AlleleComb_wi)
        
        pool.close()
        pool.join()
        
        for i in result:
            SimilarityMat[i[0]] = i[1]
        
        # Distance between two alleles, derived from similarity score

        for comb in self.AlleleComb_wo:
            # print(f"Dist: {comb}")
            distance = np.sqrt(SimilarityMat[(comb[0], comb[0])] + SimilarityMat[(comb[1], comb[1])] - 2 * SimilarityMat[comb])
            self.DistMat.loc[comb[1],comb[0]] = distance

        return

    def SaveDist(self):
        
        self.DistMat.to_csv(self.OutCSV)
        
        return