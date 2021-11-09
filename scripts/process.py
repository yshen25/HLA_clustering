#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
process pdb file downloaded from IMGT-3Dstructure database
"""

import os
from subprocess import Popen
#"B0702","B3501","B0801"
# dirs = ["B4201","B5101","B5301","B1402","B2703","B2704","B2705","B2706","B2709","B3901","B1801","B3701","B4001","B4002","B4402","B4403","B5701","B5801","B1501","B4601"]

# for InDir in dirs:
#     print(InDir)
#     for file in os.listdir(InDir):
#         if file.endswith(".pdb"):
#             os.rename(f"{InDir}/{file}", f"{InDir}/{file.split('-')[1]}")

#     os.mkdir(f"{InDir}/pdb_A")

#     for file in os.listdir(InDir):
#         if file.endswith(".pdb"):
#             extract = Popen(["extract_chain.py", f"{InDir}/{file}", "A"])
#             extract.communicate()

#     for file in os.listdir(InDir):
#         if file.endswith("_A.pdb"):
#             os.rename(f"{InDir}/{file}", f"{InDir}/pdb_A/{file}")

InDir = "B3501"
os.mkdir(f"{InDir}/pdb_A")

for file in os.listdir(InDir):
    if file.endswith(".pdb"):
        extract = Popen(["extract_chain.py", f"{InDir}/{file}", "A"])
        extract.communicate()

for file in os.listdir(InDir):
    if file.endswith("_A.pdb"):
        os.rename(f"{InDir}/{file}", f"{InDir}/pdb_A/{file}")