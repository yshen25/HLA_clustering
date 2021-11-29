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
from itertools import chain, groupby
from operator import itemgetter

import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
# from Bio.PDB import extract

from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist

from pymol import cmd

import pandas as pd

from PropertyParams import PartialCharge

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
        if atom.altloc not in [" ", "A"]:
            return 0
        else:
            return 1

def extract(structure, chain_id, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)


def PDB_trim(InDir, TemplatePDB, OutDir, OutCSV, chain="A"):
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
        print("trim:", InPDB)
        InStruct = parser.get_structure("target", f"{InDir}/{InPDB}")
        InSeq = PepBuilder.build_peptides(InStruct)[0].get_sequence()
        
        starting_loc = InStruct[0]["A"].child_list[0]._id[1]
        alignments = aligner.align(InSeq, TSeq)

        qstart = alignments[0].aligned[0][0][0] + starting_loc # alignment is 0-based, starting loc is 1-based
        #qend = alignments[0].aligned[0][-1][-1]
        qend = qstart + 178
        # use 177 to remove last amino acid of relaxed models

        # OutPDB = InPDB.split("S")[0].replace("*", "").replace(":", "_") + ".pdb"
        OutPDB = InPDB
        
        extract(InStruct, chain, qstart, qend, f"{OutDir}/{OutPDB}")
        #print(f"Trim file saved: {OutDir}/{OutPDB}, {qend-qstart+1}")
        record.append([OutPDB, qstart, qend, qend-qstart+1])

    df = pd.DataFrame(record, columns=["FileName", "qstart", "qend", "length"])
    df.to_csv(OutCSV)

    return

def alphaNbeta(InPDB):
    """
    Extract alpha helix and beta sheets aa index of input PDB (1-based)
    Return string of index range
    """
    dssp = dssp_dict_from_pdb_file(InPDB, DSSP="/home/shawn/local/dssp/mkdssp")
    secondary_structure = [dssp[0][i][1] for i in dssp[0].keys()]
    aa_index = [dssp[0][i][5] for i in dssp[0].keys()] # 1-based

    anchor_index = np.array(aa_index)[np.isin(secondary_structure, ["E", "H", "I", "G"])]

    anchor_range = []
    for _, g in groupby(enumerate(anchor_index),lambda x:x[0]-x[1]):
        group = list(map(itemgetter(1),g))
        anchor_range.append(f"{group[0]}-{group[-1]}")

    anchor_range = "+".join(anchor_range)

    return anchor_index, anchor_range

def PDB_align(InDir, refPDB, OutDir):
    """
    Superimpose query PDB to template PDB
    """
    
    cmd.load(refPDB, "template")
    _, RAchRange = alphaNbeta(refPDB)

    for InPDB in os.listdir(InDir):
        print("align:", InPDB)
        _, TAchRange = alphaNbeta(f"{InDir}/{InPDB}")

        cmd.load(f"{InDir}/{InPDB}", "target")
        cmd.align(f"target///{TAchRange}/CA", f"template///{RAchRange}/CA") # align and superimpose based on alpha helix wall and beta sheet plate

        OutPDB = f"{OutDir}/{InPDB.split('.')[0]}.pdb"
        cmd.save(OutPDB, "target")
        print(f"Align file saved: {OutPDB}")
        cmd.delete("target")
    
    return

def groove_CA_coord(PDBpath, Struct):
    """
    Input pdb and corresponding model parsed by PDBparser
    return array of coordinates of CA that form helix wall and sheet plate
    """
    anchor_index, _ = alphaNbeta(PDBpath)
    OutList = []
    for chain in Struct:
        for residue in chain:
            if residue.id[1] in anchor_index:
                atom = residue["CA"]
                X_coord, Y_coord, Z_coord = atom.coord[0:3]
                OutList.append([X_coord, Y_coord, Z_coord])
    
    return np.array(OutList)

def in_groove(PDBpath, Struct, AtomCoord, PepSurf):
    """
    input pdb, corresponding model, and atom coordinate
    return a 0/1 array of atom inside/outside of binding groove defined by helix wall and sheet plate
    """
    box = Delaunay(groove_CA_coord(PDBpath, Struct))
    InOrOut = box.find_simplex(AtomCoord)>=0
    InOrOut = InOrOut.astype(int).reshape(-1,1)

    NearPep = np.min(cdist(AtomCoord, PepSurf),1) <= 5.0
    NearPep = NearPep.astype(int).reshape(-1,1)

    return InOrOut * NearPep

def resi_depth(Struct, AtomCoord, PepSurf):
    """
    first extact the binding groove surface, then calculate the distance of each atom to binding groove
    """
    # whole protein surface
    surface = get_surface(Struct, MSMS="/home/shawn/local/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1")
    hull = Delaunay(AtomCoord)

    # remove surface vertices that fall out side of convex hull
    in_pocket = hull.find_simplex(surface)>=0
    pocket = surface[in_pocket,:]

    # remove surface vertices that is far from bound peptide
    in_groove = np.min(cdist(pocket, PepSurf),1) <= 1.0 # controls distance cut-off between peptide surface and binding groove surface
    groove = pocket[in_groove,:]

    ResiDepth = np.min(cdist(AtomCoord, groove),1).reshape(-1,1)

    return ResiDepth

def PDB_to_csv(InDir, OutDir):
    """
    Assign parameters like partial charge to each atom, and store dataframe into csv file
    """
    # with open('ffcharge.pkl', 'rb') as inf:
    #     ffcharge = pickle.load(inf)

    ffcharge = PartialCharge

    with open('pep_surf.pkl', 'rb') as inf:
        pep = pickle.load(inf)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)

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
        # i = 0 # used to record aa position
        # pep_len = len(TStruct['A'])
        for residue in TStruct['A']:
            # i += 1
            # if i >= pep_len: # if the last one
            #     break # remove the last amino acid, since ROSETTA tranforms last residue to C terminal variant
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

        OutList = np.array(OutList)
        # ResiDepth = resi_depth(TStruct, OutList[:,4:7], pep)
        # InGroove = in_groove(f"{InDir}/{InPDB}", TStruct, OutList[:,4:7], pep)
        ResiDepth = InGroove = np.zeros((OutList.shape[0],1))

        # print(OutList[:,4:7].shape)
        # print(ResiDepth.shape)

        # OutList = np.hstack((OutList, InGroove))
        # OutDF = pd.DataFrame(OutList, columns=["Residue", "ResNum", "Atom", "AtomNum", "X", "Y", "Z", "Charge", "InGroove"])

        OutList = np.hstack((OutList, ResiDepth, InGroove))
        OutDF = pd.DataFrame(OutList, columns=["Residue", "ResNum", "Atom", "AtomNum", "X", "Y", "Z", "Charge", "Depth", "InGroove"])

        OutDF.to_csv(f"{OutDir}/{InPDB.split('.')[0] + '.csv'}", index=False)
        # sys.exit()

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

def PDB_preprocess(PDBDIr, TemplatePDB, TrimDir, AlignDir, OutCSV):

    if not os.path.exists(TrimDir):
        os.makedirs(TrimDir)

    if not os.path.exists(AlignDir):
        os.makedirs(AlignDir)

    PDB_trim(PDBDIr, TemplatePDB, TrimDir, OutCSV)
    PDB_align(TrimDir, TemplatePDB, AlignDir)

    return


if __name__ == "__main__":

    # PDB_preprocess("HLAA_relax/PDB", "1i4f_Crown.pdb", "HLAA_relax/TRIM", "HLAA_relax/ALIGN", "HLAA_relax_trim.csv")
    # PDB_preprocess("HLAB_relax/PDB", "1i4f_Crown.pdb", "HLAB_relax/TRIM", "HLAB_relax/ALIGN", "HLAB_relax_trim.csv")

    # PDB_to_csv("HLAA_relax/ALIGN", "HLAA_relax/DAT")
    # PDB_to_csv("HLAB_relax/ALIGN", "HLAB_relax/DAT")
    # PDB_to_csv("crystal/A_mean/pdb", "crystal/A_mean/DAT")
    # PDB_to_csv("crystal/B_mean/pdb", "crystal/B_mean/DAT")

    # PDB_preprocess("../HLAA_pdbs", "1i4f_Crown.pdb", "../HLAA_Trimmed", "../HLAA_Aligned", "../HLAA_trim.csv")
    # PDB_preprocess("../HLAB_pdbs", "1i4f_Crown.pdb", "../HLAB_Trimmed", "../HLAB_Aligned", "../HLAB_trim.csv")
    # PDB_to_csv("../HLAA_Aligned", "../HLAA_DAT")
    # # PDB_to_csv("../HLAB_Aligned", "../HLAB_DAT")
    # CP_template("../HLAA_DAT", "../HLAA_reference_panel/DAT")
    # CP_template("../HLAB_DAT", "../HLAB_reference_panel/DAT")
    # CreateRecord("HLAA_DAT", "HLAA_rec.csv")
    # CreateRecord("HLAB_DAT", "HLAB_rec.csv")
    for allele in ["A0101", "A0201", "A3003", "A3001", "A0203", "A0205", "A0206", "A0207", "A0301", "A1101", "A6801", "A2301", "A2402"]:
    # for allele in ["B0702","B3501","B4201","B5101","B5301","B0801","B1402","B2703","B2704","B2705","B2706","B2709","B3901","B1801","B3701","B4001","B4002","B4402","B4403","B5701","B5801","B1501","B4601"]:
        PDB_preprocess(f"crystal/{allele}/pdb_A", "1i4f_Crown.pdb", f"../crystal/{allele}/TRIM", f"crystal/{allele}/ALIGN", f"{allele}_trim.csv")
    
    pass