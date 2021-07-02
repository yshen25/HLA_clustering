#!usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB import extract

from pymol import cmd
"""
parameterize the binding groove of MHC moleules to voxels
input: pdb file
output: 3D matrix
"""

def PDB_trim(InDir, TemplatePDB, OutDir):
    """
    PDB structure trim to have same length with tamplate
    """

    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=1)

    TStruct = parser.get_structure("template", TemplatePDB)
    TSeq = PepBuilder.build_peptides(TStruct)[0].get_sequence()
    
    aligner = PairwiseAligner()

    for InPDB in os.listdir(InDir):
        InStruct = parser.get_structure("target", f"{InDir}/{InPDB}")
        InSeq = PepBuilder.build_peptides(InStruct)[0].get_sequence()

        alignments = aligner.align(InSeq, TSeq)

        qstart = alignments[0].aligned[0][0][0]
        qend = alignments[0].aligned[0][-1][-1]
        #print(qstart, qend)

        OutPDB = InPDB.split(".")[0].replace("*", "").replace(":", "_") + "_trim.pdb"
        
        extract(InStruct, "A", qstart, qend, f"{OutDir}/{OutPDB}")
        print(f"Trim file saved: {OutDir}/{OutPDB}")

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

def main(PDBDIr, TemplatePDB, TrimDir, AlignDir):

    if not os.path.exists(TrimDir):
        os.makedirs(TrimDir)

    if not os.path.exists(AlignDir):
        os.makedirs(AlignDir)
    
    PDB_trim(PDBDIr, TemplatePDB, TrimDir)
    PDB_align(TrimDir, TemplatePDB, AlignDir)

    return

main("finished_pdbs", "1i4f_Crown.pdb", "Trimmed", "Aligned")