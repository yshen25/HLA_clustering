#!usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from scipy.spatial.distance import cdist
from itertools import combinations, combinations_with_replacement
from multiprocessing import Pool

import numpy as np
# from scipy.stats import hypsecant
import pandas as pd

class Calculator():

    def __init__(self, DATDir, AlleleListFile, OutCSV, depth_cut=0, ContactResi:list=None, ResiWeight:dict=None) -> None:
        
        # # fixed parameters
        # self.l = 10 # charge param
        # self.sigma = 10 # shape param
        # self.d = 0.5 # depth param
        self.OutCSV = OutCSV
        self.DATDir = DATDir

        self.ContactResi = ContactResi
        self.Weight = ResiWeight
        self.OnlyGroove = 0
        self.depth_cut = depth_cut

        with open(AlleleListFile, "r") as ALF:
            AlleleList = [line.strip() for line in ALF]
        
        # Allele combinations

        # DATList = [InDAT.split(".")[0] for InDAT in sorted(os.listdir(DATDir)) if InDAT.endswith(".csv")]
        # self.AlleleComb_wi = combinations_with_replacement(DATList, 2)
        # self.AlleleComb_wo = combinations(DATList, 2)
        self.check_files(DATDir, AlleleList)
        self.AlleleComb_wi = combinations_with_replacement(AlleleList, 2)
        self.AlleleComb_wo = combinations(AlleleList, 2)

        # Distance matrix for output, in lower triangular form

        # self.DistMat = pd.DataFrame(np.zeros((len(DATList), len(DATList))), index=DATList, columns=DATList)
        self.DistMat = pd.DataFrame(np.zeros((len(AlleleList), len(AlleleList))), index=AlleleList, columns=AlleleList)

    def check_files(self, DATDir, AlleleList):
        ExistFiles = [InDAT.split(".")[0] for InDAT in sorted(os.listdir(DATDir)) if InDAT.endswith(".csv")]
        missing = [allele for allele in AlleleList if allele not in ExistFiles]

        if missing:
            raise ValueError(f"Listed allele not exist in directory: {missing}")

        return

    def CloudSimilarity(self, CoordA, ChargeA, HydroA, CoordB, ChargeB, HydroB, WeightA, WeightB):
        """
        adjusted sum of all pair of atoms between two atom clouds
        """
        # spatial only
        # SimScore = np.sum(np.exp(-cdist(CoordA, CoordB, "euclidean")**2*0.5*self.sigma))

        # spatial + electrostatic
        # SimScore = np.sum(np.exp(-cdist(ChargeA, ChargeB, "euclidean")**2*self.l) * np.exp(-cdist(CoordA, CoordB, "euclidean")**2*0.5*self.sigma))

        # spatial + electrostatic + depth to groove
        # SimScore = np.sum(np.exp(np.reciprocal(np.outer(DepthA,DepthB))*self.d) * np.exp(-cdist(ChargeA, ChargeB, "euclidean")**2*self.l) * np.exp(-cdist(CoordA, CoordB, "euclidean")**2*0.5*self.sigma))

        # New kernel: spatial
        # SimScore = np.sum( np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))**1))

        # New kernel: spatial + electrostatic
        # SimScore = np.sum( np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))**1) * np.reciprocal(np.cosh(5*cdist(ChargeA, ChargeB, "euclidean"))**1))

        # New kernel: spatial + electrostatic + residue weight
        # SimScore = np.sum( np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))) * np.reciprocal(np.cosh(5*cdist(ChargeA, ChargeB, "euclidean"))) * np.outer(WeightA, WeightB) )

        # New kernel: spatial + electrostatic + hydrophobicity + residue weight
        SimScore = np.sum( np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))) * np.reciprocal(np.cosh(5*cdist(ChargeA, ChargeB, "euclidean"))) * np.reciprocal(np.cosh(5*cdist(HydroA, HydroB, "euclidean"))) * np.outer(WeightA, WeightB) )

        # New kernel: spatial + electrostatic + depth to groove
        # SimScore = np.sum( np.reciprocal(np.cosh(5*cdist(CoordA, CoordB, "euclidean"))**1) * np.reciprocal(np.cosh(5*cdist(ChargeA, ChargeB, "euclidean"))**1) * 1/(1 + np.exp(-(np.outer(DepthA,DepthB) - 3))) )
        # SimScore = np.sum(np.reciprocal(np.cosh(0.5*cdist(CoordA, CoordB, "euclidean"))**1) * np.reciprocal(np.cosh(5*cdist(ChargeA, ChargeB, "euclidean"))**1) * np.reciprocal(1+10*np.exp(np.outer(DepthA,DepthB) - 9)))
        return SimScore

    def ParamExtract(self, DATFile):
        DAT = pd.read_csv(DATFile)
        resnum = DAT['ResNum'].values.reshape((-1,1))
        coord = DAT[['X', 'Y', 'Z']].values
        charge = DAT['Charge'].values.reshape((-1,1))
        hydrophobicity = DAT['Hydrophobicity'].values.reshape((-1,1))
        # depth = DAT['Depth'].values.reshape((-1,1))
        # in_groove = DAT['InGroove'].values.reshape((-1,1))

        return resnum, coord, charge, hydrophobicity

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
        
        resnumA, CoordA, ChargeA, HydroA = self.ParamExtract(f"{self.DATDir}/{comb[0]}.csv")
        resnumB, CoordB, ChargeB, HydroB = self.ParamExtract(f"{self.DATDir}/{comb[1]}.csv")
        # print(CoordA.shape)

        ## ===== filter atoms =====

        # skip hydrogens
        passA = ~np.isnan(HydroA).reshape(-1)
        passB = ~np.isnan(HydroB).reshape(-1)

        resnumA = resnumA[passA]
        CoordA = CoordA[passA]
        ChargeA = ChargeA[passA]
        HydroA = HydroA[passA]

        resnumB = resnumB[passB]
        CoordB = CoordB[passB]
        ChargeB = ChargeB[passB]
        HydroB = HydroB[passB]

        # calculate only on contact residues
        if self.ContactResi:
            passA = np.isin(resnumA, self.ContactResi).reshape(-1)
            resnumA = resnumA[passA]
            CoordA = CoordA[passA]
            ChargeA = ChargeA[passA]
            HydroA = HydroA[passA]
            # DepthA = DepthA[passA]

            passB = np.isin(resnumB, self.ContactResi).reshape(-1)
            resnumB = resnumB[passB]
            CoordB = CoordB[passB]
            ChargeB = ChargeB[passB]
            HydroB = HydroB[passB]
            # DepthB = DepthB[passB]
        
        """
        if self.OnlyGroove:
            passA = GrooveA == 1
            passA = passA.reshape(-1)
            resnumA = resnumA[passA]
            CoordA = CoordA[passA]
            ChargeA = ChargeA[passA]
            # DepthA = DepthA[passA]

            passB = GrooveB == 1
            passB = passB.reshape(-1)
            resnumB = resnumB[passB]
            CoordB = CoordB[passB]
            ChargeB = ChargeB[passB]
            # DepthB = DepthB[passB]

        if self.depth_cut:
            passA = DepthA <= self.depth_cut
            passA = passA.reshape(-1)
            resnumA = resnumA[passA]
            CoordA = CoordA[passA]
            ChargeA = ChargeA[passA]
            # DepthA = DepthA[passA]

            passB = DepthB <= self.depth_cut
            passB = passB.reshape(-1)
            resnumB = resnumB[passB]
            CoordB = CoordB[passB]
            ChargeB = ChargeB[passB]
            # DepthB = DepthB[passB]
        """

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
        
        return (comb, self.CloudSimilarity(CoordA, ChargeA, HydroA, CoordB, ChargeB, HydroB, WeightA, WeightB))

    def CalcDist(self):
        
        """
        Distance of two molecules in PDB files
        """
        
        # Similarity score of two alleles
        SimilarityMat = {}
        """
        for comb in self.AlleleComb_wi:
            SimilarityMat[comb] = self.CalcSim(comb)
            SimilarityMat[comb] = pool.apply_async(self.CalcSim, comb)
        """
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

# def main(DATDir, OutCSV):
#     ContactResi = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,152,156,158,159,163,167,171]
#     # DATDir = "test"
#     # OutCSV = "test.csv"
#     calc = Calculator(DATDir, OutCSV, ContactResi=ContactResi)
#     # t1 = time.time()
#     # x = threading.Thread(target=calc.CalcDist(), args=(1,))
#     # x.start()
#     calc.CalcDist()
#     # pool = Pool(os.cpu_count())
#     # pool.map(calc.CalcDist(), 1)
#     # t2 = time.time()
#     # print(t2-t1, "seconds")

#     calc.SaveDist()
#     return

# if __name__ == '__main__':
#     main("HLAA_reference_panel/DAT", "HLAAref_distance_scd.csv")
#     main("HLAB_reference_panel/DAT", "HLABref_distance_scd.csv")