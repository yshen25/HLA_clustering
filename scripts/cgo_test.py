from pymol.cgo import ALPHA, BEGIN, END, LINE_STRIP, SPHERE, COLOR, VERTEX
from pymol import cmd

import numpy as np

# for i in range(10):
#     obj.extend([SPHERE, 0, 0, i, 1])

# obj = []

# datafiles = ["A01_01_group.csv", "A02_01_group.csv", "A03_01_group.csv", "A23_01_group.csv"]
datafiles = ["B07_02_group.csv", "B08_01_group.csv", "B14_02_group.csv", "B15_01_group.csv", "B18_01_group.csv", "B57_01_group.csv"]

for filename in datafiles:
    # with open(filename, "r") as inf:
    #     obj = inf.read().replace("\n", ",").split(",")
    #     cmd.load_cgo(obj,filename.split(".")[0])

    obj = np.genfromtxt(filename, delimiter=',')
    obj = list(obj.reshape(-1))
    cmd.load_cgo(obj,filename.split(".")[0])
# obj.extend([ALPHA, 0.25])

# obj.extend( [COLOR, 0.9, 0.9, 1.0, SPHERE, 0, 1, 0, 0.1] )
# obj.extend( [COLOR, 0.8, 0.8, 1.0, SPHERE, 0, 2, 0, 0.1] )
# obj.extend( [COLOR, 0.7, 0.7, 1.0, SPHERE, 0, 3, 0, 0.1] )
# obj.extend( [COLOR, 0.6, 0.6, 1.0, SPHERE, 0, 4, 0, 0.1] )

# obj.extend( [COLOR, 0.9, 0.9, 1.0 , SPHERE, 1, 1, 0, 0.1] )
# obj.extend( [COLOR, 0.8, 0.8, 1.0 , SPHERE, 1, 2, 0, 0.1] )
# obj.extend( [COLOR, 0.7, 0.7, 1.0 , SPHERE, 1, 3, 0, 0.1] )
# obj.extend( [COLOR, 0.6, 0.6, 1.0 , SPHERE, 1, 4, 0, 0.1] )

# obj.extend( [COLOR, 0.9, 0.9, 1.0 , SPHERE, 2, 1, 0, 0.1] )
# obj.extend( [COLOR, 0.8, 0.8, 1.0 , SPHERE, 2, 2, 0, 0.1] )
# obj.extend( [COLOR, 0.7, 0.7, 1.0 , SPHERE, 2, 3, 0, 0.1] )
# obj.extend( [COLOR, 0.6, 0.6, 1.0 , SPHERE, 2, 4, 0, 0.1] )