#!usr/bin/env python3
# -*- coding: utf-8 -*-

from re import template
from Bio import PDB
from Bio import pairwise2
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import substitution_matrices, PairwiseAligner

import time
"""
parameterize the binding groove of MHC moleules to voxels
input: pdb file
output: 3D matrix
"""

def PDB_align(InPDB, TemplatePDB):
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
    blosum62 = substitution_matrices.load("BLOSUM62")

    # s1 = time.time()
    alignments = aligner.align(InSeq, TSeq)
    # s2 = time.time()

    qstart = alignments[0].aligned[0][0][0]
    qend = alignments[0].aligned[0][-1][-1]
    print(qstart, qend)



    return

Template_PDB = "1i4f(A0201).pdb"
Input_PDB = "A0101_0078.pdb"
print("run")
PDB_align(Input_PDB, Template_PDB)