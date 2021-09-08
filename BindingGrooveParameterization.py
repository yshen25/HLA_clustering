#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
parameterize the binding groove of MHC moleules to voxels
input: pdb file
output: csv file
"""

import os
import pickle
import shutil
from string import digits

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.ResidueDepth import get_surface, min_dist

from pymol import cmd

from biopandas.pdb import PandasPdb
import pandas as pd

class ChainSelector:
    """
    Adapted from Bio.PDB.Dice module
    Only accepts residues with right chainid, between start and end.
    Remove waters and ligands. Only use model 0 by default.
    Hydrogens are kept
    """

    def __init__(self, chain_id, start, end, model_id=0):
        """Initialize the class."""
        self.chain_id = chain_id
        self.start = start
        self.end = end
        self.model_id = model_id

    def accept_model(self, model):
        """Verify if model match the model identifier."""
        # model - only keep model 0
        if model.get_id() == self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        """Verify if chain match chain identifier."""
        if chain.get_id() == self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            print(f"WARNING: Icode {icode} at position {resseq}")
        if self.start <= resseq <= self.end:
            return 1
        return 0

    def accept_atom(self, atom):
        """Modified to accept all atoms including hydrogens"""
        return 1

def extract(structure, chain_id, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)


def PDB_trim(InDir, TemplatePDB, OutDir, OutCSV):
    """
    PDB structure trim to have same length with tamplate
    """
    record = []

    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=1)

    TStruct = parser.get_structure("template", TemplatePDB)
    TSeq = PepBuilder.build_peptides(TStruct)[0].get_sequence()
    
    aligner = PairwiseAligner()
    aligner.gap_score = -10 # no gap wanted

    for InPDB in os.listdir(InDir):
        InStruct = parser.get_structure("target", f"{InDir}/{InPDB}")
        InSeq = PepBuilder.build_peptides(InStruct)[0].get_sequence()

        alignments = aligner.align(InSeq, TSeq)

        qstart = alignments[0].aligned[0][0][0]
        #qend = alignments[0].aligned[0][-1][-1]
        qend = qstart + 178

        OutPDB = InPDB.split(".")[0].replace("*", "").replace(":", "_") + "_trim.pdb"
        
        extract(InStruct, "A", qstart, qend, f"{OutDir}/{OutPDB}")
        #print(f"Trim file saved: {OutDir}/{OutPDB}, {qend-qstart+1}")
        record.append([OutPDB, qstart, qend, qend-qstart+1])

    df = pd.DataFrame(record, columns=["FileName", "qstart", "qend", "length"])
    df.to_csv(OutCSV)

    return

def PDB_align(InDir, TemplatePDB, OutDir):
    """
    Superimpose query PDB to template PDB
    """
    
    cmd.load(TemplatePDB, "template")

    for InPDB in os.listdir(InDir):
        cmd.load(f"{InDir}/{InPDB}", "target")
        cmd.align("target", "template")

        OutPDB = f"{OutDir}/{InPDB.split('.')[0]}_align.pdb"
        cmd.save(OutPDB, "target")
        print(f"Align file saved: {OutPDB}")
        cmd.delete("target")
    
    return

def PDB_to_csv(InDir, OutDir):
    import sys
    """
    Assign parameters like partial charge to each atom, and store dataframe into pickle file
    """
    with open('ffcharge.pkl', 'rb') as inf:
        ffcharge = pickle.load(inf)

    for InPDB in os.listdir(InDir):
        print(InPDB)
        OutList = []
        '''
        OutPKL = InPDB.split(".")[0] + ".pdb"
        Struct_df.to_pickle(OutPKL)
        print(f"Pickle file saved: {OutDir}/{OutPKL}")
        '''
        parser = PDBParser(PERMISSIVE=1)

        TStruct = parser.get_structure("struct", f"{InDir}/{InPDB}")[0]
        # surface = get_surface(TStruct, MSMS="/home/shawn/work_bench/HLA_clustering/MSMS/msms.x86_64Linux2.2.6.1.staticgcc") # protein surface
        # print(surface)
        # from matplotlib import pyplot
        # from mpl_toolkits.mplot3d import Axes3D

        """fig = pyplot.figure()
        ax = Axes3D(fig)
        ax.scatter(surface[:,0], surface[:,1], surface[:,2])
        pyplot.show()
        sys.exit()"""
        for chain in TStruct:
            for residue in chain:
                #print(residue.__dict__)
                #print(dir(residue))
                ResName = residue.resname
                ResNum = residue.id[1]
                
                if ResName == "HIS": # distinguish three protonation states of histidine
                    
                    if residue.has_id("HE2"):
                        
                        if residue.has_id("HD1"):
                            ResName = "HIP" # H on both epsilon and delta N
                        else:
                            ResName = "HIE" # H on epsilon N only

                    else:
                        ResName = "HID" # H on delta N only. By default HIS is HID
                
                for atom in residue:
                    #print(atom.__dict__)
                    X_coord, Y_coord, Z_coord = atom.coord[0:3]
                    OutList.append([ResName, ResNum, atom.name, atom.serial_number, X_coord, Y_coord, Z_coord
                    , ffcharge[ResName][atom.name.lstrip(digits)]])
                    #, min_dist(atom.coord, surface)

        OutDF = pd.DataFrame(OutList, columns=["Residue", "ResNum", "Atom", "AtomNum", "X", "Y", "Z", "Charge"])
        OutDF.to_csv(f"{OutDir}/{InPDB.split('S')[0] + '.csv'}", index=False)
        #sys.exit()

    return

def CP_template(DATDir, RefDir):
    RefList = ["A01_01","A02_01","A02_02","A02_03","A02_04","A02_05","A02_06","A02_07","A02_14","A02_17"
    ,"A03_01","A11_01","A23_01","A24_02","A26_01","A26_02","A26_03","A29_02","A30_01","A30_02","A30_03"
    ,"A30_04","A31_01","A32_01","A33_01","A33_03","A66_01","A68_01","A68_02","A69_01","A74_01","B07_02"
    ,"B07_03","B07_05","B08_01","B08_02","B14_02","B15_01","B15_02","B15_03","B15_08","B15_09","B15_10"
    ,"B15_12","B15_13","B15_16","B15_17","B15_18","B18_01","B27_02","B27_03","B27_04","B27_05","B27_06"
    ,"B27_07","B27_09","B35_01","B35_03","B37_01","B38_01","B39_01","B39_02","B39_09","B40_01","B40_02"
    ,"B40_06","B42_01","B44_02","B44_03","B45_01","B46_01","B48_01","B51_01","B51_02","B51_03","B52_01"
    ,"B53_01","B54_01","B55_01","B55_02","B56_01","B57_01","B57_02","B58_01","B58_02","B67_01","B73_01"
    ,"B78_01"]

    for InDAT in os.listdir(DATDir):
        if InDAT.split(".")[0] in RefList:
            shutil.copy2(f"{DATDir}/{InDAT}", RefDir)
            print(f"Copy reference allele: {DATDir}/{InDAT}")

    return

def CreateRecord(DATDir, RecFile):
    
    RefList = ["A01_01","A02_01","A02_02","A02_03","A02_04","A02_05","A02_06","A02_07","A02_14","A02_17"
    ,"A03_01","A11_01","A23_01","A24_02","A26_01","A26_02","A26_03","A29_02","A30_01","A30_02","A30_03"
    ,"A30_04","A31_01","A32_01","A33_01","A33_03","A66_01","A68_01","A68_02","A69_01","A74_01","B07_02"
    ,"B07_03","B07_05","B08_01","B08_02","B14_02","B15_01","B15_02","B15_03","B15_08","B15_09","B15_10"
    ,"B15_12","B15_13","B15_16","B15_17","B15_18","B18_01","B27_02","B27_03","B27_04","B27_05","B27_06"
    ,"B27_07","B27_09","B35_01","B35_03","B37_01","B38_01","B39_01","B39_02","B39_09","B40_01","B40_02"
    ,"B40_06","B42_01","B44_02","B44_03","B45_01","B46_01","B48_01","B51_01","B51_02","B51_03","B52_01"
    ,"B53_01","B54_01","B55_01","B55_02","B56_01","B57_01","B57_02","B58_01","B58_02","B67_01","B73_01"
    ,"B78_01"]

    record = []
    for InDAT in sorted(os.listdir(DATDir)):
        if InDAT.split(".")[0] in RefList:
            record.append([InDAT.split(".")[0], 1])
        else:
            record.append([InDAT.split(".")[0], 0])
    df = pd.DataFrame(record, columns=["Allele", "reference"])
    df.to_csv(RecFile, index=False)

    return

def main(PDBDIr, TemplatePDB, TrimDir, AlignDir, OutCSV):

    if not os.path.exists(TrimDir):
        os.makedirs(TrimDir)

    if not os.path.exists(AlignDir):
        os.makedirs(AlignDir)

    PDB_trim(PDBDIr, TemplatePDB, TrimDir, OutCSV)
    PDB_align(TrimDir, TemplatePDB, AlignDir)

    return

# main("HLAA_pdbs", "1i4f_Crown.pdb", "HLAA_Trimmed", "HLAA_Aligned", "HLAA_trim.csv")
# main("HLAB_pdbs", "1i4f_Crown.pdb", "HLAB_Trimmed", "HLAB_Aligned", "HLAB_trim.csv")
PDB_to_csv("HLAA_Aligned", "HLAA_DAT")
PDB_to_csv("HLAB_Aligned", "HLAB_DAT")
# CP_template("HLAA_DAT", "HLAA_reference_panel")
# CP_template("HLAB_DAT", "HLAB_reference_panel")
# CreateRecord("HLAA_DAT", "HLAA_rec.csv")
# CreateRecord("HLAB_DAT", "HLAB_rec.csv")