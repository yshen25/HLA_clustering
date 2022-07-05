#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract atom and residue coordinates from protein structures in pdb format
The atomic or coarse-grained residue coordinates are stored in csv files
"""
import os
import glob
from string import digits
from itertools import combinations
from pathlib import Path
from typing import Union
from itertools import groupby
from operator import itemgetter
import shutil

import numpy as np
# from scipy.stats import hypsecant
import pandas as pd
from pymol import cmd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import matplotlib.cm as cm
import open3d as o3d

from .PropertyParams import AtomicMass, PartialCharge, AtomicHydrophobicity

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
        if chain.get_id() in self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, loc, icode = residue.get_id()
        chain_ = residue.parent.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            print(f"WARNING: Icode {icode} at position {loc}")
        if self.start[self.chain_id.index(chain_)] <= loc <= self.end[self.chain_id.index(chain_)]:
            return 1
        return 0

    def accept_atom(self, atom):
        """Modified to accept all atoms including hydrogens"""
        
        # _hydrogen = re.compile("[123 ]*H.*")

        # if _hydrogen.match(atom.get_id()): # remove hydrogens
        #     return 0
        if "H" in atom.get_id(): # new way to remove hydrogens
            return 0

        if atom.altloc not in [" ", "A"]: # remove altloc atoms
            return 0

        return 1

def extract(structure, chain_id, start, end, filename):
    """Write out selected portion of structure to <filename (pdb)>."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)

    return

def PDB_renumber(struct, start:list):
    """
    renumber pdb files, the first residue in pdb file is indexed as start number
    each chain is given a start number in alphabat order
    """
    for i, chain in enumerate(struct[0]):
    #loop 1: renumber residues to negative number to avoid errors
        residue_id = -1
        for residue in chain.get_residues():
            residue.id=(' ',residue_id,' ')
            residue_id -= 1
        #loop 2
        residue_id = start[i]
        for residue in chain.get_residues():
            #print(chain.get_id(), residue_id)
            residue.id=(' ',residue_id,' ')
            residue_id += 1
    return struct

def PDB_trim(InDir, TemplatePDB, OutDir, chain="A", length=[179], template_start_id=[2]):
    """
    PDB structure trim to have same residues as the tamplate
    """
    # length=[179] # HLA1
    # start_id=[2] # HLA1, 1i4f
    # length = [80, 85] # HLA2
    # template_start_id = [4,7] # HLA2, 1klu
    # record = []

    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    TStruct = parser.get_structure("template", TemplatePDB)
    TSeqs = PepBuilder.build_peptides(TStruct)
    
    aligner = PairwiseAligner()
    aligner.gap_score = -10 # no gap wanted

    for InPDB in os.listdir(InDir):
        if InPDB.endswith(".pdb"):
            print("trim:", InPDB)
            InStruct = parser.get_structure("target", f"{InDir}/{InPDB}")
            Seqs = PepBuilder.build_peptides(InStruct[0])
            qends = []
            qstarts = []
            for i, chain_identifier in enumerate(chain):
                InSeq = Seqs[i].get_sequence()
                TSeq = TSeqs[i].get_sequence()
            
                starting_loc = InStruct[0][chain_identifier].child_list[0]._id[1] # index of the first residue
                alignments = aligner.align(InSeq, TSeq)

            ## === archive, stand-alone trim ===
            # qstart = alignments[0].aligned[0][0][0] + starting_loc # alignment is 0-based, starting loc is 1-based
            # qend = qstart + 178 # fixed length
            #qend = alignments[0].aligned[0][-1][-1] # aligned portion
            
            ## === use with PDB_renumber ===
                # qstart = starting_loc - alignments[0].aligned[0][0][0] + template_start_id[i] -1
                qstart = starting_loc + alignments[0].aligned[1][0][0] - alignments[0].aligned[0][0][0] + template_start_id[i] -1
                # starting loc is calibarated, for the 1st residue in template is not loc 1
                qend = length[i]

                qstarts.append(qstart)
                qends.append(qend)
                
                # record.append([InPDB, chain_identifier, qstart, qend, qend-qstart+1])

            OutPDB = InPDB

            # InStruct = PDB_renumber(InStruct, qstart)
            # extract(InStruct, chain, start_id, qends, f"{OutDir}/{OutPDB}") # HLA1
            InStruct = PDB_renumber(InStruct, qstarts)
            extract(InStruct, chain, [2,2], qends, f"{OutDir}/{OutPDB}")
            #print(f"Trim file saved: {OutDir}/{OutPDB}, {qend-qstart+1}")
            

    # df = pd.DataFrame(record, columns=["FileName", "chain", "qstart", "qend", "length"])
    # df.to_csv(OutCSV)

    return

def PDB_trim(InDir, TemplatePDB, OutDir):
    """
    PDB structure trim to have same residues as the tamplate
    """
    
    return

def BB_RMSD(Qpdb, Rpdb):

    cmd.load(Qpdb, "query")
    cmd.load(Rpdb, "ref")

    rmsd = cmd.align(f"query////CA", f"ref////CA", cycles=0, transform=0)[0] # backbone RMSD

    cmd.delete("query")
    cmd.delete("ref")

    return rmsd

def centroid_structure(InPDBDir):
    """
    find the structure that has least RMSD with other structures
    used as alignment reference
    """
    # list all pdb files in input dir
    InFiles = glob.glob(f"{InPDBDir}/*.pdb")
    
    # initialize RMSD matrix
    RMSD_Mat = pd.DataFrame(np.zeros((len(InFiles), len(InFiles)), dtype=float), index=InFiles, columns=InFiles)

    # calculate pairwise RMSD and fill in the matrix
    comb = combinations(InFiles, 2)
    for pair in comb:
        RMSD_Mat.loc[pair[1],pair[0]] = BB_RMSD(pair[1], pair[0])

    RMSD_Mat = RMSD_Mat.add(RMSD_Mat.T, fill_value=0) # lower triangle to square
    # print(RMSD_Mat)
    RMSD_sum = RMSD_Mat.sum(axis=1)
    return RMSD_sum.idxmin()

def alphaNbeta(InPDB, dssp_path=None):
    """
    Extract alpha helix and beta sheets aa index of input PDB (1-based)
    Return string of index range
    ====================
    InPDB: Input PDB file
    dssp_path: Path to mkdssp excutable
    --------------------
    output
    anchor_index: Indexes of residues that forms alpha helix and beta strands
    anchor_range: Ranges of indexes
    """
    if dssp_path:
        if not os.path.exists(dssp_path):
            raise ValueError("DSSP not found!")
    else:
        dssp_path = shutil.which("mkdssp")
        if not dssp_path:
            raise ValueError("DSSP not found!")
    
    dssp = dssp_dict_from_pdb_file(InPDB, DSSP=dssp_path)
    secondary_structure = [dssp[0][i][1] for i in dssp[0].keys()]
    aa_index = [dssp[0][i][5] for i in dssp[0].keys()] # 1-based

    anchor_index = np.array(aa_index)[np.isin(secondary_structure, ["E", "H", "I", "G"])]

    anchor_range = []
    for _, g in groupby(enumerate(anchor_index),lambda x:x[0]-x[1]):
        group = list(map(itemgetter(1),g))
        anchor_range.append(f"{group[0]}-{group[-1]}")

    anchor_range = "+".join(anchor_range)

    return anchor_index, anchor_range

def PDB_align(InDir, OutDir, refPDB=None, SSAlign=False):
    """
    Superimpose input PDBs to the reference PDB
    if not specified, the centroid structure is used
    """
    # Find the centroid structure as reference
    if refPDB is None:
        print("Reference structure is not specified\n\tSearching for centroid structure...")
        refPDB = centroid_structure(InDir)
        print(f"Centroid: {refPDB}")

    if SSAlign:
        _, RAchRange = alphaNbeta(refPDB)

    cmd.load(refPDB, "template")

    for InPDB in glob.glob(f"{InDir}/*.pdb"):
        # print("align:", InPDB)
        cmd.load(InPDB, "target")
        cmd.h_add(selection="(resn 'GLY' and name 'CA')") # add hydrogen atoms for Glycine, especifically for crystal structures
        cmd.alter("(resn 'GLY' and name 'H01')", "name='1HA'")
        cmd.alter("(resn 'GLY' and name 'H02')", "name='2HA'")
        cmd.alter("(resn 'GLY' and name 'H03')", "name='2HA'")

        if SSAlign: # align and superimpose based on secondary structures
            _, TAchRange = alphaNbeta(InPDB)
            cmd.align(f"target///{TAchRange}/CA", f"template///{RAchRange}/CA")
        
        else:
            cmd.align(f"target////CA", f"template////CA") # align and superimpose

        cmd.save(f"{OutDir}/{Path(InPDB).name}", "target")
        # print(f"Align file saved: {OutPDB}")
        cmd.delete("target")
    
    cmd.delete("template")
    return

def PDB_to_csv(InDir, OutDir):
    """
    Assign parameters like partial charge to each atom, and store dataframe into csv file
    """
    # with open('ffcharge.pkl', 'rb') as inf:
    #     ffcharge = pickle.load(inf)

    ffcharge = PartialCharge
    ffhydro = AtomicHydrophobicity

    # with open('pep_surf.pkl', 'rb') as inf:
    #     pep = pickle.load(inf)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)

    for InPDB in glob.glob(f"{InDir}/*.pdb"):
        OutList = []

        parser = PDBParser(QUIET=True)

        TStruct = parser.get_structure("struct", InPDB)[0]
        # i = 0 # used to record aa position
        # pep_len = len(TStruct['A'])
        for chain in TStruct:
            for residue in chain:
                # i += 1
                # if i >= pep_len: # if the last one
                #     break # remove the last amino acid, since ROSETTA tranforms last residue to C terminal variant
                ResName = residue.resname
                ResNum = residue.id[1]
                
                if ResName == "HIS": # distinguish three protonation states of histidine
                    
                    if residue.has_id("HE2"):
                        
                        if residue.has_id("HD1"):
                            ResName2 = "HIP" # H on both epsilon and delta N
                        else:
                            ResName2 = "HIE" # H on epsilon N only

                    else:
                        ResName2 = "HID" # H on delta N only. By default HIS is HID

                elif ResName == "MSE": # no partial charge value for MSE
                    ResName = ResName2 = "MET"
                
                else:
                    ResName2 = ResName
                
                for atom in residue:
                    # Rosetta relax will automaticly change last Oxigen from O to OXT
                    if atom.name == "OXT":
                        atom.name = "O"
                    # change Se in MSE to S
                    elif atom.name == "SE":
                        atom.name = "SD"

                    X_coord, Y_coord, Z_coord = atom.coord[0:3]
                    OutList.append([ResName, chain.id, ResNum, atom.name, atom.serial_number, X_coord, Y_coord, Z_coord
                    , ffcharge[ResName2][atom.name.lstrip(digits)], ffhydro[ResName][atom.name.lstrip(digits)], 1])

        OutList = np.array(OutList)

        OutDF = pd.DataFrame(OutList, columns=["Residue", "Chain", "ResNum", "Atom", "AtomNum", "X", "Y", "Z", "Charge", "Hydrophobicity", "Accept"])

        OutDF.to_csv(f"{OutDir}/{Path(InPDB).stem}.csv", index=False)

    return

def CoarseGrain(InDAT, OutCGDAT):
    """
    center-of-mass of side chain
    """
    # if OutFile is None:
    #     OutFile = DATFile.split(".")[0] + "_CG.csv"

    FullAtom_df = pd.read_csv(InDAT)
    # DAT = DAT[~DAT.Atom.isin(["N", "CA", "C", "O"])] # filter out backbone atoms

    ResGen = FullAtom_df.groupby(by=['Chain', 'ResNum'])

    CG_row = []
    for resi in ResGen:
        x = y = z = 0
        total_mass = 0
        chain = resi[0][0]
        resnum = resi[0][1]
        resatoms_df = resi[1]
        # print(resatoms)
        resname = resatoms_df["Residue"].iloc[0]

        for atom in resatoms_df[["Atom", "X", "Y", "Z"]].to_numpy():
            # print(atom)
            if atom[0] in ["N", "CA", "C", "O"]:
                continue # filter out backbone atoms

            if not atom[0][0].isdigit():
                mass = AtomicMass[atom[0][0]] # real atom name is the first letter of pdb atom name

            else:
                mass = AtomicMass[atom[0][1]]

            total_mass += mass
            x += mass * atom[1]
            y += mass * atom[2]
            z += mass * atom[3]

        if total_mass == 0:
            CG_row.append([chain, resnum, resname, np.nan, np.nan, np.nan, 0])
        else:
            CG_row.append([chain, resnum, resname, x/total_mass, y/total_mass, z/total_mass, 1])

    CG_DAT = pd.DataFrame(CG_row, columns=["Chain", "ResNum", "Residue", "X", "Y", "Z", "Weight"])
    CG_DAT.to_csv(OutCGDAT, index=False)

    return

def PointCloudCG(DATDir, OutDir):
    """
    Input: full atom DAT directory
    Output: coarse-grained DAT file
    """
    for InDAT in glob.glob(f"{DATDir}/*.csv"):
        CoarseGrain(InDAT, f"{OutDir}/{Path(InDAT).stem}_CG.csv")

    return

def element_depth():
    return

def depth_threshold():
    return

def min_dist():
    return

def dist_threshold():
    return

def PDBtoPointCloud(InputDir:Union[str, bytes, os.PathLike], AlignDir:Union[str, bytes, os.PathLike], 
    PCDir:Union[str, bytes, os.PathLike], TrimDir:Union[str, bytes, os.PathLike]=None, CGDir:Union[str, bytes, os.PathLike]=None, 
    SSAlign:bool=False, RefPDB:Union[str, bytes, os.PathLike]=None)->None:
    """
    Turn PDB file into atom cloud
    ====================
    Input: Input directory that contains PDB files
    SSAlign: Do structure alignment based on residues that form secondary structures, alphahelixes and beta strands
    trim: Trim every input PDB files according to Reference PDB, the Reference must be specified
    RefPDB: PDB file used as reference in structure alignment. If not specified, the centroid structure will be calculated and used
    --------------------
    Output: Trimmed and aligned PDB files, CSV pointcloud files
    """
    # ===== validate parameters =====
    # If trim, the RefPDB must be specified
    # if trim, the TRIM dir must be specified
    if TrimDir:
        if not RefPDB:
            raise ValueError("If enable trim, a RefPDB must be specified")
        if not os.path.exists(TrimDir):
            print(f"Create directory for trimmed PDBs: {TrimDir}")
        PDB_trim()
    
    # ===== Align =====
    if not os.path.exists(AlignDir):
        print(f"Create directory for aligned PDBs: {AlignDir}")
        os.mkdir(AlignDir)
    PDB_align(InputDir, AlignDir, refPDB=RefPDB, SSAlign=SSAlign)

    # ===== To point cloud =====
    if not os.path.exists(PCDir):
        print(f"Create directory for point cloud CSVs: {PCDir}")
        os.mkdir(PCDir)
    PDB_to_csv(AlignDir, PCDir)

    # ===== Coarse graining =====
    if CGDir:
        if not os.path.exists(CGDir):
            print(f"Create directory for coarse-grained point cloud CSVs: {CGDir}")
            os.mkdir(CGDir)
        PointCloudCG(PCDir, CGDir)
    return

class assign_weight:
    def read_CG_DAT(self, InCG_DAT):
        self.filename = InCG_DAT
        df = pd.read_csv(InCG_DAT)
        self.df = df.set_index(['Chain', 'ResNum'])
        return

    def by_resi_number(self, weight_dict):
        reform = [[i,j,weight_dict[i][j]] for i in weight_dict.keys() for j in weight_dict[i].keys()]
        self.weight = pd.DataFrame(reform, columns=["Chain", "ResNum", "Weight"])
        self.weight = self.weight.set_index(['Chain', 'ResNum'])
        
        return
    def update_weight(self):
        # initialize by giving 0 in every position
        self.df['Weight'] = 0
        
        # residue weight as multi-indexed dataframe
        # self.weight_dict_to_df(weight_dict)

        # change weight via update
        self.df.update(self.weight)
        return

    def save_CG_DAT(self):
        self.df.to_csv(self.filename)
        return

def Show_PointCloud(DATDir):
    pointcloudlist = []
    cmap = cm.ScalarMappable(cmap="coolwarm")
    for InDAT in glob.glob(f"{DATDir}/*.csv"):
        df = pd.read_csv(InDAT)
        P1 = o3d.geometry.PointCloud()
        coord = df[['X','Y','Z']].to_numpy()
        weight = df['Weight'].to_numpy()

        color1 = cmap.to_rgba(weight)[:,0:3]
        P1.points = o3d.utility.Vector3dVector(coord)
        P1.colors = o3d.utility.Vector3dVector(color1)
        pointcloudlist.append(P1)

    o3d.visualization.draw_geometries(pointcloudlist)
    return

def CGDAT_dir_reweight(CGDAT_dir, weight_dict):
    re_weight = assign_weight()
    re_weight.by_resi_number(weight_dict)

    re_weight.weight
    for CGDAT_file in glob.glob(f"{CGDAT_dir}/*.csv"):

        re_weight.read_CG_DAT(CGDAT_file)
        re_weight.update_weight()
        re_weight.save_CG_DAT()