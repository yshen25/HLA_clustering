#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download crystal structure of A*0201 and extract the comprehensive surface of bound peptide
"""
import sys
import urllib.request
import subprocess
import os
import pickle

import numpy as np
import pandas as pd

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.ResidueDepth import get_surface

import open3d as o3d

import alphashape

from BindingGrooveParameterization import PDB_align

# from pymol import cmd

# Download pdb file from IMGT-3Dstructures database

pdb_id = ["3mre","3d25","6trn","5n1y","1i4f","6o53","4u6y","5ddh","6o4z","2gtw","6o51"
    ,"6o4y","2v2w","2v2x","2vll","3bgm","3gso","3pwn","3v5h","3bh8","3fqn","3kla","3pwl"
    ,"3utq","3qfd","4u6x","2git","2gtz","3bh9","3fqr","3fqx","3o3d","3pwj","2gt9","6opd"
    ,"1duz","1i7u","1tvb","1tvh","3fqt","3fqu","3ft2","3gsu","3o3a","3gsw","1jf1","3o3e"
    ,"5enw","3h7b","5mer","3myj","1t1z","2guo","2x4q","3ft4","3gjf","3gsv","3o3b","3rew"
    ,"5hhp","6pte","3fqw","2av1","3ft3","3gsr","1t1y","2clr","2x70","3giv","3hpj","3v5d"
    ,"4k7f","5swq","5hhn","2av7","2x4u","3gsx","3ixa","5hhq","3gsq","1duy","1jht","6ptb"
    ,"1t21","1i1y","1i7r","1im3","1qew","1s8d","1s9w","1t1w","1t1x","1t20","1t22","3bhb"
    ,"3i6g","3mgt","6apn","1eey","1eez","1s9y","2x4o","2x4p","2x4r","2x4t","3mgo","3v5k"
    ,"2x4n","2c7u","1qr1","1b0g","1hhi","1hhj","1hhk","1s9x","5hhm","2x4s","1hhg","3hla"
    ,"3to2","1akj","5mep","1i1f","1i7t","3i6k","1b0r","3gjg","3hae","5hho","1hhh","4wuu"
    ,"6nca","1p7q","1hla","3h9h","3mrg","3mrb","3mrk","3mrr","3mrd","3mrc","3mrj","3mrm"
    ,"3mr9","3mri","3mrp","3mrq","3mrf","3mrn","3mro","3mrh","3mrl","4jfo","4jfp","4jfq"
    ,"2j8u","2uwe","2jcc","2p5e","2p5w","3hg1","1lp9","3utt","3uts","2bnq","2bnr","2f54"
    ,"2f53","3qfj","2gj6","1ao7","3pwp","3h9s","1qrn","1qse","1qsf","3d3v","3d39","6rsy"
    ,"3o4l","1bd2","2pye","3qeq","3qdm","3qdj","3qdg","6tro","1oga","2vlr","2vlj","2vlk"
        ,"6rpb","6rpa","6rp9","3gsn","5tez"]

PDBDIr = "peptide_pdb"
TrimDir = "peptide_trim"
AlignDir = "peptide_aln"
PepDir = "peptide_pep"
TemplatePDB = "1i4f_Crown.pdb"
OutCSV = "peptide_trim.csv"

if not os.path.exists(TrimDir):
    os.makedirs(TrimDir)

if not os.path.exists(AlignDir):
    os.makedirs(AlignDir)

if not os.path.exists(PepDir):
    os.makedirs(PepDir)

def step1():
    i=1
    for id in pdb_id:
        print(f"Downloading: {id}  {i}/{len(pdb_id)}")
        urllib.request.urlretrieve(f"http://www.imgt.org/3Dstructure-DB/IMGT-FILE/IMGT-{id.upper()}.pdb.gz", f"{PDBDIr}/{id}.pdb.gz")
        i += 1
    
    return

# decompress donwloaded files
def step2():
    files = [f"{PDBDIr}/{id}.pdb.gz" for id in pdb_id]
    process = subprocess.Popen(["gunzip"]+files)
    process.communicate()

    return

# extract A (alpha subunit) and C (peptide) chain, and trim A chain
class ChainSelector:
    """
    Adapted from Bio.PDB.Dice module
    Only accepts residues with right chainid, between start and end.
    Remove waters and ligands. Only use model 0 by default.
    Hydrogens are kept
    """

    def __init__(self, start, end, model_id=0):
        """Initialize the class."""
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
        if chain.get_id() in ("A", "C"):
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # # residue - between start and end
        hetatm_flag, loc, icode = residue.get_id()
        chain = residue.get_parent().get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            # print(f"WARNING: Icode {icode} at position {loc}")
            return 0
        # print(resseq)
        if chain == "A" and (self.start <= loc <= self.end):
            return 1
        elif chain == "C":
            return 1

        return 0
        # return 1

    def accept_atom(self, atom):
        """Modified to accept all atoms including hydrogens"""
        return 1

def extract(structure, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)

def renumber(struct):
    #loop 1: renumber residues to negative number to avoid errors
    chain_id = ""
    for residue in struct.get_residues():
        hetatm_flag, _, icode = residue.get_id()
        if icode != ' ':
            
            continue # skip alternative positions

        chain = residue.get_parent()
        if chain_id != chain.get_id():
            chain_id = chain.get_id()
            residue_id = -1
        #print(chain.get_id(), residue_id)
        residue.id=(hetatm_flag, residue_id, icode)
        residue_id -= 1
    #loop 2: reverse the negative number
    for residue in struct.get_residues():
        hetatm_flag, _, icode = residue.get_id()
        residue.id = (hetatm_flag, -residue.get_id()[1], icode)

    return struct

def step3():

    record = []

    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=1, QUIET=1)

    TStruct = parser.get_structure("template", TemplatePDB)
    TSeq = PepBuilder.build_peptides(TStruct)[0].get_sequence()

    aligner = PairwiseAligner()
    aligner.gap_score = -10 # no gap wanted

    for InPDB in os.listdir(PDBDIr):

        InStruct = parser.get_structure("target", f"{PDBDIr}/{InPDB}")

        InStruct = renumber(InStruct)
        
        '''
        try:
            InPep = PepBuilder.build_peptides(InStruct[0]['A'])
            InPep2 = PepBuilder.build_peptides(InStruct[0]['C'])
            print(InPDB, len(InPep[0].get_sequence()), len(InPep2[0].get_sequence()))

        except:
            print("!!!", InPDB)
        '''

        InPep = PepBuilder.build_peptides(InStruct[0]['A'])
        InSeq = InPep[0].get_sequence()

        alignments = aligner.align(InSeq, TSeq)

        qstart = alignments[0].aligned[0][0][0]
        # qend = alignments[0].aligned[0][-1][-1]
        qend = qstart + 178

        OutPDB = InPDB.split(".")[0] + "_trim.pdb"
        
        extract(InStruct, qstart, qend, f"{TrimDir}/{OutPDB}")

        record.append([OutPDB, qstart, qend, qend-qstart+1])
        # sys.exit()

    df = pd.DataFrame(record, columns=["FileName", "qstart", "qend", "length"])
    df.to_csv(OutCSV)

    return

### align to A0201 template, then remove A chain
def step4():

    PDB_align(TrimDir, "1i4f_Crown.pdb", AlignDir)
    parser = PDBParser(PERMISSIVE=1, QUIET=1)

    for InPDB in os.listdir(AlignDir):
        print(InPDB)
        # Inparser = PDBParser(PERMISSIVE=1, QUIET=0)
        InStruct = parser.get_structure("target", f"{AlignDir}/{InPDB}")

        OutPDB = InPDB.split("_")[0] + "_C.pdb"
        
        extract(InStruct, -1, -1, f"{PepDir}/{OutPDB}")

    return

### calculate the surface of all peptides
### then get the alpha-shape of all surface vertex points as comprehensive surface

def step5():

    pep_surf_bag = np.empty((0,3), float)
    parser = PDBParser(PERMISSIVE=1, QUIET=1)
    for InPDB in os.listdir(PepDir):
        # print(InPDB)
        InStruct = parser.get_structure("target", f"{PepDir}/{InPDB}")
        surface = get_surface(InStruct[0], MSMS="/home/shawn/local/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1")
        # print(surface, type(surface))
        pep_surf_bag = np.vstack((pep_surf_bag, surface))
        # sys.exit()
    
    print(pep_surf_bag.shape)

    alpha_shape = alphashape.alphashape(pep_surf_bag, alpha=0.9).vertices

    # visualize peptide surfaces

    # P1 = o3d.geometry.PointCloud()
    # P1.points = o3d.utility.Vector3dVector(pep_surf_bag)
    # P1.paint_uniform_color([1, 0, 0])

    # P2 = o3d.geometry.PointCloud()
    # P2.points = o3d.utility.Vector3dVector(alpha_shape)
    # P2.paint_uniform_color([0, 0, 1])

    # o3d.visualization.draw_geometries([P1, P2])

    # save peptide surface in pickle

    with open("pep_surf.pkl", "wb") as ouf:
        pickle.dump(alpha_shape, ouf)
    
    return

step5()