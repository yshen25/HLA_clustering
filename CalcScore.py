#!usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
from scipy.spatial.distance import cdist

from biopandas.pdb import PandasPdb
import pandas as pd

def CloudSimilarity(MatA, MatB, sigma):
    """
    adjusted sum of all pair of atoms between two atom clouds
    """
    SumDist = np.sum(np.exp(-cdist(MatA, MatB)/2*sigma*sigma))
    return SumDist

def CalDist(TargetPDB, ReferencePDB, sigma):
    """
    Distance of two molecules in PDB files
    """

    TCoord = PandasPdb().read_pdb(TargetPDB).df["ATOM"][["x_coord", "y_coord", "z_coord"]]
    RCoord = PandasPdb().read_pdb(ReferencePDB).df["ATOM"][["x_coord", "y_coord", "z_coord"]]

    Dist = np.sqrt(CloudSimilarity(TCoord, TCoord, sigma)
                   +CloudSimilarity(RCoord, RCoord, sigma)
                   -2*CloudSimilarity(TCoord, RCoord, sigma))

    return Dist

#print(CalDist("Aligned/A01_03S_0099_trim_align.pdb", "1i4f_Crown.pdb", 0.1))

def main():
    RefPanel = os.listdir("reference_panel")
    TargetPDBs = os.listdir("Aligned")
    ScoreList = []
    for target in TargetPDBs:
        AlleleScore = []
        for Ref in RefPanel:
            AlleleScore.append(CalDist(f"Aligned/{target}", f"reference_panel/{Ref}", 0.1))

        ScoreList.append(AlleleScore)

    df = pd.DataFrame(ScoreList, columns=RefPanel)
    df.to_csv("distance.csv")
    return

main()