#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
module for using and parsing msms
"""
from subprocess import Popen, PIPE, DEVNULL
import tempfile
import numpy as np

from scipy.spatial.distance import cdist

class SolventAccessible():
    """
    predict if provided points are inside or outside of the solvent accessible volume of a protein
    """
    def __init__(self, MSMS_Dir=None, pdb_to_xyzr=None, MSMS=None) -> None:
        
        # path to msms package
        if  MSMS_Dir:
            self.MSMS_Dir = MSMS_Dir
        else:
            self.MSMS_Dir = "/Users/ys0/local/msms_x86_64Darwin_2.6.1"
        
        if pdb_to_xyzr:
            self.pdb_to_xyzr = pdb_to_xyzr
        else:
            self.pdb_to_xyzr = "pdb_to_xyzr"
        
        if MSMS:
            self.MSMS = MSMS
        else:
            self.MSMS = "msms.x86_64Darwin.2.6.1"

    def run_msms(self, InPDB):
        """
        input: PDB file name
        output: msms output filename including reduced surface points and normals
        """
        pdb_to_xyzr_run = Popen([f"{self.MSMS_Dir}/{self.pdb_to_xyzr}", InPDB], stdout=PIPE, stderr=PIPE)
        xyzr_out, error = pdb_to_xyzr_run.communicate()

        if error:
            raise ValueError(error)
        
        xyzr_temp = tempfile.NamedTemporaryFile()

        with open(xyzr_temp.name, "w") as fh:
            fh.write(xyzr_out.decode())

        OutFile = tempfile.NamedTemporaryFile()

        msms_run = Popen([f"{self.MSMS_Dir}/{self.MSMS}", "-if", xyzr_temp.name, "-of", OutFile.name], stdout=DEVNULL)
        msms_run.communicate()
        
        return self.parse_msms(OutFile.name)

    def parse_msms(self, msmsFile):
        """
        input: msms output file name (.vert file)
        output: surface points and normals
        """
        with open(msmsFile) as fp:
            vertex_list = []
            norm_list = []
            for l in fp:
                sl = l.split()
                if len(sl) != 9:
                    # skip header
                    continue
                vl = [float(x) for x in sl[0:3]]
                nl = [float(x) for x in sl[3:6]]
                vertex_list.append(vl)
                norm_list.append(nl)

        vertex_list = np.array(vertex_list)
        norm_list = np.array(norm_list)

        return vertex_list, norm_list

    def InorOut(self, GridList:np.ndarray, TargetPDB:str) -> np.ndarray:
        """
        input: numpy array of 3D points
        output: numpy array of points inside
        algorithm: calculate dot product of vector from nearest surface point to 
            query point and normal of nearest surface point, if > 0 is outside
        """
        vertexL, normL = self.run_msms(TargetPDB)
        nearest_point_index = np.argmin(cdist(GridList, vertexL), axis=1)
        
        nearest_point = vertexL[nearest_point_index]
        nearest_normal = normL[nearest_point_index]

        projection = np.dot(GridList - nearest_point, nearest_normal)

        outside_bool = projection > 0

        return GridList[outside_bool]