#!usr/bin/env python3
# -*- coding: utf-8 -*-

from re import template
from Bio import PDB
from Bio import pairwise2
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB import extract

from pymol import cmd

import time
"""
parameterize the binding groove of MHC moleules to voxels
input: pdb file
output: 3D matrix
"""

def PDB_trim(InPDB, TemplatePDB):
    """
    PDB structure pre-processing
    """
    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=1)
    InStruct = parser.get_structure("target", InPDB)
    TStruct = parser.get_structure("template", TemplatePDB)

    InSeq = PepBuilder.build_peptides(InStruct)[0].get_sequence()
    TSeq = PepBuilder.build_peptides(TStruct)[0].get_sequence()
    
    aligner = PairwiseAligner()

    # s1 = time.time()
    alignments = aligner.align(InSeq, TSeq)
    # s2 = time.time()

    qstart = alignments[0].aligned[0][0][0]
    qend = alignments[0].aligned[0][-1][-1]
    print(qstart, qend)

    OutPDB = ".".split(InPDB)[0] + "_trim.pdb"
    extract(InStruct, "A", start, stop, OutPDB)
    print(f"Trim file saved: {OutPDB}")

    return OutPDB

def PDB_align(InPDB, TemplatePDB):
    InName = ".".split(InPDB)[0]
    TName = ".".split(TemplatePDB)[0]
    cmd.load(InPDB, InName)
    cmd.load(TemplatePDB, TName)
    cmd.align(InName, TName)

    OutPDB = f"{InName}_on_{TName}.pdb"
    cmd.save(OutPDB, InName)
    print(f"Align file saved: {OutPDB}")
    return OutPDB


Template_PDB = "1i4f(A0201).pdb"
Input_PDB = "A0101_0078.pdb"
print("run")
trim_file = PDB_trim(Input_PDB, Template_PDB)
align_file = PDB_align(trim_file, Template_PDB)