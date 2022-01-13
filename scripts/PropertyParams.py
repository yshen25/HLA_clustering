#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines atomic properties like mass and partial charge
For other scripts to import and lookup
"""
import pandas as pd
import numpy as np

AtomicMass = {'H':1.00794, 'C':12.011, 'N':14.00674, 'O':15.9994, 'S':32.066}
# Based on 1993 IUPAC table of standard atomic weights of the elements

PartialCharge = {
    "ALA":{
        "N":-0.415700,
        "H":0.271900,
        "CA":0.033700,
        "HA":0.082300,
        "CB":-0.182500,
        "HB":0.060300,
        "C":0.597300,
        "O":-0.567900
        },
    "ARG":{
        "N":-0.347900,
        "H":0.274700,
        "CA":-0.263700,
        "HA":0.156000,
        "CB":-0.000700,
        "HB":0.032700,
        "CG":0.039000,
        "HG":0.028500,
        "CD":0.048600,
        "HD":0.068700,
        "NE":-0.529500,
        "HE":0.345600,
        "CZ":0.807600,
        "NH1":-0.862700,
        "NH2":-0.862700,
        "HH1":0.447800,
        "HH2":0.447800,
        "C":0.734100,
        "O":-0.589400
        },
    "ASN":{
        "N":-0.415700,
        "H":0.271900,
        "CA":0.014300,
        "HA":0.104800,
        "CB":-0.204100,
        "HB":0.079700,
        "CG":0.713000,
        "OD1":-0.593100,
        "ND2":-0.919100,
        "HD2":0.419600,
        "C":0.597300,
        "O":-0.567900
        },
    "ASP":{
        "N":-0.516300,
        "H":0.293600,
        "CA":0.038100,
        "HA":0.088000,
        "CB":-0.030300,
        "HB":-0.012200,
        "CG":0.799400,
        "OD1":-0.801400,
        "OD2":-0.801400,
        "C":0.536600,
        "O":-0.581900
        },
    "CYS":{
        "N":-0.415700,
        "H":0.271900,
        "CA":0.021300,
        "HA":0.112400,
        "CB":-0.123100,
        "HB":0.111200,
        "SG":-0.311900,
        "HG":0.193300,
        "C":0.597300,
        "O":-0.567900
        },
    "GLN":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.003100,
        "HA":0.085000,
        "CB":-0.003600,
        "HB":0.017100,
        "CG":-0.064500,
        "HG":0.035200,
        "CD":0.695100,
        "OE1":-0.608600,
        "NE2":-0.940700,
        "HE2":0.425100,
        "C":0.597300,
        "O":-0.567900
        },
    "GLU":{
        "N":-0.516300,
        "H":0.293600,
        "CA":0.039700,
        "HA":0.110500,
        "CB":0.056000,
        "HB":-0.017300,
        "CG":0.013600,
        "HG":-0.042500,
        "CD":0.805400,
        "OE1":-0.818800,
        "OE2":-0.818800,
        "C":0.536600,
        "O":-0.581900
        },
    "GLY":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.025200,
        "HA":0.069800,
        "C":0.597300,
        "O":-0.567900
        },
    "HID":{
        "N":-0.415700,
        "H":0.271900,
        "CA":0.018800,
        "HA":0.088100,
        "CB":-0.046200,
        "HB":0.040200,
        "CG":-0.026600,
        "ND1":-0.381100,
        "HD1":0.364900,
        "HD2":0.114700,
        "CE1":0.205700,
        "HE1":0.139200,
        "NE2":-0.572700,
        "CD2":0.129200,
        "C":0.597300,
        "O":-0.567900
        },
    "HIE":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.058100,
        "HA":0.136000,
        "CB":-0.007400,
        "HB":0.036700,
        "CG":0.186800,
        "ND1":-0.543200,
        "CE1":0.163500,
        "HE1":0.143500,
        "NE2":-0.279500,
        "HE2":0.333900,
        "CD2":-0.220700,
        "HD2":0.186200,
        "C":0.597300,
        "O":-0.567900
        },
    "HIP":{
        "N":-0.347900,
        "H":0.274700,
        "CA":-0.135400,
        "HA":0.121200,
        "CB":-0.041400,
        "HB2":0.081000,
        "HB3":0.081000,
        "CG":-0.001200,
        "ND1":-0.151300,
        "HD1":0.386600,
        "CE1":-0.017000,
        "HE1":0.268100,
        "NE2":-0.171800,
        "HE2":0.391100,
        "CD2":-0.114100,
        "HD2":0.231700,
        "C":0.734100,
        "O":-0.589400
        },
    "ILE":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.059700,
        "HA":0.086900,
        "CB":0.130300,
        "HB":0.018700,
        "CG1":-0.043000,
        "CG2":-0.320400,
        "HG1":0.023600,
        "HG2":0.088200,
        "CD1":-0.066000,
        "HD1":0.018600,
        "C":0.597300,
        "O":-0.567900
        },
    "LEU":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.051800,
        "HA":0.092200,
        "CB":-0.110200,
        "HB":0.045700,
        "CG":0.353100,
        "HG":-0.036100,
        "CD1":-0.412100,
        "CD2":-0.412100,
        "HD1":0.100000,
        "HD2":0.100000,
        "C":0.597300,
        "O":-0.567900
        },
    "LYS":{
        "N":-0.347900,
        "H":0.274700,
        "CA":-0.240000,
        "HA":0.142600,
        "CB":-0.009400,
        "HB":0.036200,
        "CG":0.018700,
        "HG":0.010300,
        "CD":-0.047900,
        "HD":0.062100,
        "CE":-0.014300,
        "HE":0.113500,
        "NZ":-0.385400,
        "HZ":0.340000,
        "C":0.734100,
        "O":-0.589400
        },
    "MET":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.023700,
        "HA":0.088000,
        "CB":0.034200,
        "HB":0.024100,
        "CG":0.001800,
        "HG":0.044000,
        "SD":-0.273700,
        "CE":-0.053600,
        "HE":0.068400,
        "C":0.597300,
        "O":-0.567900
        },
    "PHE":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.002400,
        "HA":0.097800,
        "CB":-0.034300,
        "HB":0.029500,
        "CG":0.011800,
        "CD1":-0.125600,
        "CD2":-0.125600,
        "HD1":0.133000,
        "HD2":0.133000,
        "CE1":-0.170400,
        "CE2":-0.170400,
        "HE1":0.143000,
        "HE2":0.143000,
        "CZ":-0.107200,
        "HZ":0.129700,
        "C":0.597300,
        "O":-0.567900
        },
    "PRO":{
        "N":-0.254800,
        "CD":0.019200,
        "HD":0.039100,
        "CG":0.018900,
        "HG":0.021300,
        "CB":-0.007000,
        "HB":0.025300,
        "CA":-0.026600,
        "HA":0.064100,
        "C":0.589600,
        "O":-0.574800
        },
    "SER":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.024900,
        "HA":0.084300,
        "CB":0.211700,
        "HB":0.035200,
        "OG":-0.654600,
        "HG":0.427500,
        "C":0.597300,
        "O":-0.567900
        },
    "THR":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.038900,
        "HA":0.100700,
        "CB":0.365400,
        "HB":0.004300,
        "CG2":-0.243800,
        "HG1":0.410200,
        "HG2":0.064200,
        "OG1":-0.676100,
        "C":0.597300,
        "O":-0.567900
        },
    "TRP":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.027500,
        "HA":0.112300,
        "CB":-0.005000,
        "HB":0.033900,
        "CG":-0.141500,
        "CD1":-0.163800,
        "CD2":0.124300,
        "HD1":0.206200,
        "NE1":-0.341800,
        "HE1":0.341200,
        "HE3":0.170000,
        "CE2":0.138000,
        "CE3":-0.238700,
        "CZ2":-0.260100,
        "CZ3":-0.197200,
        "CH2":-0.113400,
        "HH2":0.141700,
        "HZ2":0.157200,
        "HZ3":0.144700,
        "C":0.597300,
        "O":-0.567900
        },
    "TYR":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.001400,
        "HA":0.087600,
        "CB":-0.015200,
        "HB":0.029500,
        "CG":-0.001100,
        "CD1":-0.190600,
        "CD2":-0.190600,
        "CE1":-0.234100,
        "CE2":-0.234100,
        "HE1":0.165600,
        "HE2":0.165600,
        "CZ":0.322600,
        "OH":-0.557900,
        "HH":0.399200,
        "HD1":0.169900,
        "HD2":0.169900,
        "C":0.597300,
        "O":-0.567900
        },
    "VAL":{
        "N":-0.415700,
        "H":0.271900,
        "CA":-0.087500,
        "HA":0.096900,
        "CB":0.298500,
        "HB":-0.029700,
        "CG1":-0.319200,
        "CG2":-0.319200,
        "HG1":0.079100,
        "HG2":0.079100,
        "C":0.597300,
        "O":-0.567900
        }
}
# Derived from amino12.lib in Amber forcefield, which is used by ff14SB method

# AASim = {
#     "SER":{"SER":1.00, "ARG":0.49, "LEU":0.33, "PRO":0.66, "THR":0.73, "ALA":0.54, "VAL":0.42, "GLY":0.74, "ILE":0.34, "PHE":0.28, "TYR":0.33, "CYS":0.48, "HIS":0.59, "GLN":0.68, "ASN":0.79, "LYS":0.44, "ASP":0.70, "GLU":0.63, "MET":0.37, "TRP":0.18},
#     "ARG":{"SER":0.49, "ARG":1.00, "LEU":0.53, "PRO":0.52, "THR":0.67, "ALA":0.48, "VAL":0.55, "GLY":0.42, "ILE":0.55, "PHE":0.55, "TYR":0.64, "CYS":0.16, "HIS":0.87, "GLN":0.80, "ASN":0.60, "LYS":0.88, "ASP":0.55, "GLU":0.75, "MET":0.58, "TRP":0.53},
#     "LEU":{"SER":0.33, "ARG":0.53, "LEU":1.00, "PRO":0.54, "THR":0.57, "ALA":0.55, "VAL":0.85, "GLY":0.36, "ILE":0.98, "PHE":0.90, "TYR":0.83, "CYS":0.08, "HIS":0.54, "GLN":0.47, "ASN":0.29, "LYS":0.50, "ASP":0.20, "GLU":0.36, "MET":0.93, "TRP":0.72},
#     "PRO":{"SER":0.66, "ARG":0.52, "LEU":0.54, "PRO":1.00, "THR":0.82, "ALA":0.87, "VAL":0.68, "GLY":0.80, "ILE":0.56, "PHE":0.47, "TYR":0.49, "CYS":0.21, "HIS":0.64, "GLN":0.65, "ASN":0.58, "LYS":0.52, "ASP":0.50, "GLU":0.57, "MET":0.60, "TRP":0.32},
#     "THR":{"SER":0.73, "ARG":0.67, "LEU":0.57, "PRO":0.82, "THR":1.00, "ALA":0.73, "VAL":0.68, "GLY":0.73, "ILE":0.59, "PHE":0.52, "TYR":0.57, "CYS":0.31, "HIS":0.78, "GLN":0.80, "ASN":0.70, "LYS":0.64, "ASP":0.60, "GLU":0.70, "MET":0.62, "TRP":0.40},
#     "ALA":{"SER":0.54, "ARG":0.48, "LEU":0.55, "PRO":0.87, "THR":0.73, "ALA":1.00, "VAL":0.70, "GLY":0.72, "ILE":0.56, "PHE":0.47, "TYR":0.48, "CYS":0.09, "HIS":0.60, "GLN":0.58, "ASN":0.48, "LYS":0.51, "ASP":0.41, "GLU":0.50, "MET":0.61, "TRP":0.31},
#     "VAL":{"SER":0.42, "ARG":0.55, "LEU":0.85, "PRO":0.68, "THR":0.68, "ALA":0.70, "VAL":1.00, "GLY":0.49, "ILE":0.87, "PHE":0.77, "TYR":0.74, "CYS":0.11, "HIS":0.61, "GLN":0.55, "ASN":0.38, "LYS":0.55, "ASP":0.29, "GLU":0.44, "MET":0.90, "TRP":0.59},
#     "GLY":{"SER":0.74, "ARG":0.42, "LEU":0.36, "PRO":0.80, "THR":0.73, "ALA":0.72, "VAL":0.49, "GLY":1.00, "ILE":0.37, "PHE":0.29, "TYR":0.32, "CYS":0.26, "HIS":0.54, "GLN":0.60, "ASN":0.63, "LYS":0.41, "ASP":0.56, "GLU":0.54, "MET":0.41, "TRP":0.14},
#     "ILE":{"SER":0.34, "ARG":0.55, "LEU":0.98, "PRO":0.56, "THR":0.59, "ALA":0.56, "VAL":0.87, "GLY":0.37, "ILE":1.00, "PHE":0.90, "TYR":0.85, "CYS":0.08, "HIS":0.56, "GLN":0.49, "ASN":0.31, "LYS":0.53, "ASP":0.22, "GLU":0.38, "MET":0.95, "TRP":0.72},
#     "PHE":{"SER":0.28, "ARG":0.55, "LEU":0.90, "PRO":0.47, "THR":0.52, "ALA":0.47, "VAL":0.77, "GLY":0.29, "ILE":0.90, "PHE":1.00, "TYR":0.90, "CYS":0.05, "HIS":0.53, "GLN":0.46, "ASN":0.27, "LYS":0.53, "ASP":0.18, "GLU":0.35, "MET":0.87, "TRP":0.81},
#     "TYR":{"SER":0.33, "ARG":0.64, "LEU":0.83, "PRO":0.49, "THR":0.57, "ALA":0.48, "VAL":0.74, "GLY":0.32, "ILE":0.85, "PHE":0.90, "TYR":1.00, "CYS":0.10, "HIS":0.61, "GLN":0.54, "ASN":0.33, "LYS":0.60, "ASP":0.26, "GLU":0.43, "MET":0.83, "TRP":0.83},
#     "CYS":{"SER":0.48, "ARG":0.16, "LEU":0.08, "PRO":0.21, "THR":0.31, "ALA":0.09, "VAL":0.11, "GLY":0.26, "ILE":0.08, "PHE":0.05, "TYR":0.10, "CYS":1.00, "HIS":0.19, "GLN":0.28, "ASN":0.35, "LYS":0.06, "ASP":0.28, "GLU":0.21, "MET":0.09, "TRP":0.00},
#     "HIS":{"SER":0.59, "ARG":0.87, "LEU":0.54, "PRO":0.64, "THR":0.78, "ALA":0.60, "VAL":0.61, "GLY":0.54, "ILE":0.56, "PHE":0.53, "TYR":0.61, "CYS":0.19, "HIS":1.00, "GLN":0.89, "ASN":0.68, "LYS":0.85, "ASP":0.62, "GLU":0.81, "MET":0.60, "TRP":0.47},
#     "GLN":{"SER":0.68, "ARG":0.80, "LEU":0.47, "PRO":0.65, "THR":0.80, "ALA":0.58, "VAL":0.55, "GLY":0.60, "ILE":0.49, "PHE":0.46, "TYR":0.54, "CYS":0.28, "HIS":0.89, "GLN":1.00, "ASN":0.79, "LYS":0.75, "ASP":0.72, "GLU":0.87, "MET":0.53, "TRP":0.40},
#     "ASN":{"SER":0.79, "ARG":0.60, "LEU":0.29, "PRO":0.58, "THR":0.70, "ALA":0.48, "VAL":0.38, "GLY":0.63, "ILE":0.31, "PHE":0.27, "TYR":0.33, "CYS":0.35, "HIS":0.68, "GLN":0.79, "ASN":1.00, "LYS":0.56, "ASP":0.89, "GLU":0.80, "MET":0.34, "TRP":0.19},
#     "LYS":{"SER":0.44, "ARG":0.88, "LEU":0.50, "PRO":0.52, "THR":0.64, "ALA":0.51, "VAL":0.55, "GLY":0.41, "ILE":0.53, "PHE":0.53, "TYR":0.60, "CYS":0.06, "HIS":0.85, "GLN":0.75, "ASN":0.56, "LYS":1.00, "ASP":0.53, "GLU":0.74, "MET":0.56, "TRP":0.49},
#     "ASP":{"SER":0.70, "ARG":0.55, "LEU":0.20, "PRO":0.50, "THR":0.60, "ALA":0.41, "VAL":0.29, "GLY":0.56, "ILE":0.22, "PHE":0.18, "TYR":0.26, "CYS":0.28, "HIS":0.62, "GLN":0.72, "ASN":0.89, "LYS":0.53, "ASP":1.00, "GLU":0.79, "MET":0.26, "TRP":0.16},
#     "GLU":{"SER":0.63, "ARG":0.75, "LEU":0.36, "PRO":0.57, "THR":0.70, "ALA":0.50, "VAL":0.44, "GLY":0.54, "ILE":0.38, "PHE":0.35, "TYR":0.43, "CYS":0.21, "HIS":0.81, "GLN":0.87, "ASN":0.80, "LYS":0.74, "ASP":0.79, "GLU":1.00, "MET":0.41, "TRP":0.29},
#     "MET":{"SER":0.37, "ARG":0.58, "LEU":0.93, "PRO":0.60, "THR":0.62, "ALA":0.61, "VAL":0.90, "GLY":0.41, "ILE":0.95, "PHE":0.87, "TYR":0.83, "CYS":0.09, "HIS":0.60, "GLN":0.53, "ASN":0.34, "LYS":0.56, "ASP":0.26, "GLU":0.41, "MET":1.00, "TRP":0.69},
#     "TRP":{"SER":0.18, "ARG":0.53, "LEU":0.72, "PRO":0.32, "THR":0.40, "ALA":0.31, "VAL":0.59, "GLY":0.14, "ILE":0.72, "PHE":0.81, "TYR":0.83, "CYS":0.00, "HIS":0.47, "GLN":0.40, "ASN":0.19, "LYS":0.49, "ASP":0.16, "GLU":0.29, "MET":0.69, "TRP":1.00}
# }
AASim_matrix = [
    [1,0.49,0.33,0.66,0.73,0.54,0.42,0.74,0.34,0.28,0.33,0.48,0.59,0.68,0.79,0.44,0.7,0.63,0.37,0.18],
    [0.49,1,0.53,0.52,0.67,0.48,0.55,0.42,0.55,0.55,0.64,0.16,0.87,0.8,0.6,0.88,0.55,0.75,0.58,0.53],
    [0.33,0.53,1,0.54,0.57,0.55,0.85,0.36,0.98,0.9,0.83,0.08,0.54,0.47,0.29,0.5,0.2,0.36,0.93,0.72],
    [0.66,0.52,0.54,1,0.82,0.87,0.68,0.8,0.56,0.47,0.49,0.21,0.64,0.65,0.58,0.52,0.5,0.57,0.6,0.32],
    [0.73,0.67,0.57,0.82,1,0.73,0.68,0.73,0.59,0.52,0.57,0.31,0.78,0.8,0.7,0.64,0.6,0.7,0.62,0.4],
    [0.54,0.48,0.55,0.87,0.73,1,0.7,0.72,0.56,0.47,0.48,0.09,0.6,0.58,0.48,0.51,0.41,0.5,0.61,0.31],
    [0.42,0.55,0.85,0.68,0.68,0.7,1,0.49,0.87,0.77,0.74,0.11,0.61,0.55,0.38,0.55,0.29,0.44,0.9,0.59],
    [0.74,0.42,0.36,0.8,0.73,0.72,0.49,1,0.37,0.29,0.32,0.26,0.54,0.6,0.63,0.41,0.56,0.54,0.41,0.14],
    [0.34,0.55,0.98,0.56,0.59,0.56,0.87,0.37,1,0.9,0.85,0.08,0.56,0.49,0.31,0.53,0.22,0.38,0.95,0.72],
    [0.28,0.55,0.9,0.47,0.52,0.47,0.77,0.29,0.9,1,0.9,0.05,0.53,0.46,0.27,0.53,0.18,0.35,0.87,0.81],
    [0.33,0.64,0.83,0.49,0.57,0.48,0.74,0.32,0.85,0.9,1,0.1,0.61,0.54,0.33,0.6,0.26,0.43,0.83,0.83],
    [0.48,0.16,0.08,0.21,0.31,0.09,0.11,0.26,0.08,0.05,0.1,1,0.19,0.28,0.35,0.06,0.28,0.21,0.09,0],
    [0.59,0.87,0.54,0.64,0.78,0.6,0.61,0.54,0.56,0.53,0.61,0.19,1,0.89,0.68,0.85,0.62,0.81,0.6,0.47],
    [0.68,0.8,0.47,0.65,0.8,0.58,0.55,0.6,0.49,0.46,0.54,0.28,0.89,1,0.79,0.75,0.72,0.87,0.53,0.4],
    [0.79,0.6,0.29,0.58,0.7,0.48,0.38,0.63,0.31,0.27,0.33,0.35,0.68,0.79,1,0.56,0.89,0.8,0.34,0.19],
    [0.44,0.88,0.5,0.52,0.64,0.51,0.55,0.41,0.53,0.53,0.6,0.06,0.85,0.75,0.56,1,0.53,0.74,0.56,0.49],
    [0.7,0.55,0.2,0.5,0.6,0.41,0.29,0.56,0.22,0.18,0.26,0.28,0.62,0.72,0.89,0.53,1,0.79,0.26,0.16],
    [0.63,0.75,0.36,0.57,0.7,0.5,0.44,0.54,0.38,0.35,0.43,0.21,0.81,0.87,0.8,0.74,0.79,1,0.41,0.29],
    [0.37,0.58,0.93,0.6,0.62,0.61,0.9,0.41,0.95,0.87,0.83,0.09,0.6,0.53,0.34,0.56,0.26,0.41,1,0.69],
    [0.18,0.53,0.72,0.32,0.4,0.31,0.59,0.14,0.72,0.81,0.83,0,0.47,0.4,0.19,0.49,0.16,0.29,0.69,1],
]
AA_order = ["SER","ARG","LEU","PRO","THR","ALA","VAL","GLY","ILE","PHE","TYR","CYS","HIS","GLN","ASN","LYS","ASP","GLU","MET","TRP"]
AASim = pd.DataFrame(AASim_matrix, columns=AA_order, index=AA_order)
# 1-D/215, D is Grantham's distance

# Nicolau DV Jr. et al. 2014
# 12 atom types, DGwif
AtomicHydrophobicity = {
    "ALA":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169
        },
    "GLY":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "HA":np.nan
        },
    "SER":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":-0.1217,
        "OG":0.4881
        },
    "THR":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":-0.1217,
        "OG1":0.4881,
        "CG2":0.2169
        },
    "LEU":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":1.5299,
        "CD1":-0.9227,
        "CD2":-0.9227
        },
    "ILE":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG1":0.2169,
        "CG2":-0.9227,
        "CD1":0.2169
        },
    "VAL":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":1.5299,
        "CG1":-0.9227,
        "CG2":-0.9227
        },
    "ASN":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.2645,
        "OD1":0.4881,
        "ND2":-0.0062
        },
    "GLN":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":0.2169,
        "CD":-0.2645,
        "OE1":0.4881,
        "NE2":-0.0062
        },
    "ARG":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":0.2169,
        "CD":-0.1217,
        "NE":-0.0062,
        "CZ":-0.2645,
        "NH1":0.3544,
        "NH2":0.3544
        },
    "HIS":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.1217,
        "ND1":0.3544,
        "CD2":-0.1217,
        "CE1":-0.1217,
        "NE2":0.3544
        },
    "TRP":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.2607,
        "CD1":-0.1217,
        "CD2":-0.2607,
        "NE1":-0.0062,
        "CE2":-0.1217,
        "CE3":-0.2607,
        "CZ2":-0.2607,
        "CZ3":-0.2607,
        "CH2":-0.2607
        },
    "PHE":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.2607,
        "CD1":-0.2607,
        "CD2":-0.2607,
        "CE1":-0.2607,
        "CE2":-0.2607,
        "CZ":-0.2607
        },
    "TYR":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.2607,
        "CD1":-0.2607,
        "CD2":-0.2607,
        "CE1":-0.2607,
        "CE2":-0.2607,
        "CZ":-0.1217,
        "OH":0.4881
        },
    "GLU":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":0.2169,
        "CD":-0.2645,
        "OE1":0.7653,
        "OE2":0.7653
        },
    "ASP":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.2645,
        "OD1":0.7653,
        "OD2":0.7653
        },
    "LYS":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":0.2169,
        "CD":0.2169,
        "CE":-0.1217,
        "NZ":-0.1231
        },
    "PRO":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":0.2169,
        "CD":-0.1217
        },
    "CYS":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":-0.1217,
        "SG":0.4989
        },
    "MET":{
        "N":-0.0062,
        "CA":-0.1217,
        "C":-0.2645,
        "O":0.4881,
        "CB":0.2169,
        "CG":-0.1217,
        "SD":0.4989,
        "CE":-0.1217
        }
}