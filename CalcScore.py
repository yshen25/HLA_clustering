#!usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
from scipy.spatial.distance import cdist
from itertools import combinations, combinations_with_replacement
from multiprocessing import Pool

import numpy as np

import pandas as pd

class Calculator():

    def __init__(self, DATDir, OutCSV) -> None:
        
        # fixed parameters
        self.l = 0.1 # charge param
        self.sigma = 0.1 # shape param
        self.OutCSV = OutCSV
        self.DATDir = DATDir
        
        # Allele combinations
        DATList = [a.split(".")[0] for a in sorted(os.listdir(DATDir))]
        self.AlleleComb_wi = combinations_with_replacement(DATList, 2)
        self.AlleleComb_wo = combinations(DATList, 2)
        
        # Distance matrix for output, in lower triangular form
        self.DistMat = pd.DataFrame(np.zeros((len(DATList), len(DATList))), index=DATList, columns=DATList)

    def CloudSimilarity(self, CoordA, ChargeA, CoordB, ChargeB):
        """
        adjusted sum of all pair of atoms between two atom clouds
        """
        # spatial only
        SimScore = np.sum(np.exp(-cdist(CoordA, CoordB)/2*self.sigma*self.sigma))

        # spatial + electrostatic
        # SimScore = np.sum(np.exp(-cdist(ChargeA, ChargeB, "cityblock")**2/self.l) * np.exp(-cdist(CoordA, CoordB)/2*self.sigma*self.sigma))

        return SimScore

    def ParamExtract(self, DATFile):
        DAT = pd.read_csv(DATFile)
        coord = DAT[['X', 'Y', 'Z']].values
        charge = DAT['Charge'].values.reshape((-1,1))
        return coord, charge

    def CalcSim(self, comb):
        print(f"Simi: {comb}")
        
        CoordA, ChargeA = self.ParamExtract(f"{self.DATDir}/{comb[0]}.csv")
        CoordB, ChargeB = self.ParamExtract(f"{self.DATDir}/{comb[1]}.csv")
        
        return (comb, self.CloudSimilarity(CoordA, ChargeA, CoordB, ChargeB))

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

        result = pool.map(self.CalcSim, self.AlleleComb_wi)
        
        pool.close()
        pool.join()
        
        for i in result:
            SimilarityMat[i[0]] = i[1]
        
        # Distance between two alleles, derived from similarity score

        for comb in self.AlleleComb_wo:
            print(f"Dist: {comb}")
            distance = np.sqrt(SimilarityMat[(comb[0], comb[0])] + SimilarityMat[(comb[1], comb[1])] - 2 * SimilarityMat[comb])
            self.DistMat.loc[comb[1],comb[0]] = distance

        return

    def SaveDist(self):
        
        self.DistMat.to_csv(self.OutCSV)
        
        return

def main():
    DATDir = "HLAB_DAT"
    OutCSV = "HLAB_distance_s.csv"
    # DATDir = "test"
    # OutCSV = "test.csv"
    calc = Calculator(DATDir, OutCSV)
    # t1 = time.time()
    # x = threading.Thread(target=calc.CalcDist(), args=(1,))
    # x.start()
    calc.CalcDist()
    # pool = Pool(os.cpu_count())
    # pool.map(calc.CalcDist(), 1)
    # t2 = time.time()
    # print(t2-t1, "seconds")

    calc.SaveDist()
    return

main()