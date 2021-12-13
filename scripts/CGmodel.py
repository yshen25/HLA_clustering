#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From full-atom DAT file to coarse-grained DAT file
Residues are represented by side-chain center-of-mass
"""
import os
import sys
import time
from typing_extensions import ParamSpec
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
        resatoms_df = resi[1]
        # print(resatoms)
        resname = resatoms_df["Residue"].iloc[0]

        for atom in resatoms_df[["Atom", "X", "Y", "Z"]].to_numpy():
            print(atom)
            mass = AtomicMass[atom[0][0]] # real atom name is the first letter of pdb atom name
            total_mass += mass
            x += mass * atom[1]
            y += mass * atom[2]
            z += mass * atom[3]


        CG_list.append([resnum, resname, x/total_mass, y/total_mass, z/total_mass])

    CG_DAT = pd.DataFrame(CG_list, columns=["ResNum", "Residue", "X", "Y", "Z"])
    CG_DAT.to_csv(OutFile, index=False)

    return

class CGCalculator():
    """
    Calculate distance between CG models
    """
    def __init__(self, DATDir, OutCSV, ContactResi:list=None, ResiWeight:dict=None, Pairwise:bool=False) -> None:
        
        # first, check if pairwise mode seeting is correct
        # number of residues between two molecules must be the same
        if Pairwise:
            if not ContactResi:
                raise ValueError("Must specify contact residues to enable pairwise mode!!")
        self.pairwise = Pairwise # controls pairwise matrices or all-to-all matrices
        
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
        """
        Residue similarity
        """
        if self.pairwise:
            ResiPairComb = list(zip(ResiA, ResiB))
        else:
            ResiPairComb = [(x, y) for x in ResiA for y in ResiB]
        
        ResiPairSim = np.array([AASim[i][j] for i,j in ResiPairComb])

        # print(ResiPairComb.shape)
        if self.pairwise:
            
            return ResiPairSim

        else:

            return ResiPairSim.reshape((len(ResiA), len(ResiB)))


    def CloudSimilarity(self, CoordA, ResiA, CoordB, ResiB, WeightA, WeightB):
        
        ResiPairSim_score = self.ResiPairSim(ResiA, ResiB)

        if self.pairwise:
            SimScore = np.sum(np.linalg.norm(CoordA - CoordB, ord=2, axis=1) * ResiPairSim_score * np.multiply(WeightA, WeightB))

        else:
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
        
        # print(len(CoordA), len(CoordB))
        if len(resnumA) != len(resnumB):
            print(f"{comb[0]} != {comb[1]}")
        sys.exit()
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