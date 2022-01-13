from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd

"""
draw a rectangular cuboid represeting sampling region
"""

X_range = (-16, 15)
Y_range = (-7, 7)
Z_range = (0, 10)

vertice1 = [X_range[0], Y_range[0], Z_range[0]]
vertice2 = [X_range[0], Y_range[0], Z_range[1]]
vertice3 = [X_range[0], Y_range[1], Z_range[0]]
vertice4 = [X_range[0], Y_range[1], Z_range[1]]
vertice5 = [X_range[1], Y_range[0], Z_range[0]]
vertice6 = [X_range[1], Y_range[0], Z_range[1]]
vertice7 = [X_range[1], Y_range[1], Z_range[0]]
vertice8 = [X_range[1], Y_range[1], Z_range[1]]

# in pymol.cgo file, CYLINDER = 9.0
# pymol cgo object format: [shape, starting_point(xyz), ending_point(xyz), radius, starting_color(rgb), ending_color(rgb)]

obj = [
   CYLINDER, vertice1[0], vertice1[1], vertice1[2], vertice2[0], vertice2[1], vertice2[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice1[0], vertice1[1], vertice1[2], vertice3[0], vertice3[1], vertice3[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice1[0], vertice1[1], vertice1[2], vertice5[0], vertice5[1], vertice5[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice2[0], vertice2[1], vertice2[2], vertice4[0], vertice4[1], vertice4[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice2[0], vertice2[1], vertice2[2], vertice6[0], vertice6[1], vertice6[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice3[0], vertice3[1], vertice3[2], vertice4[0], vertice4[1], vertice4[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice3[0], vertice3[1], vertice3[2], vertice7[0], vertice7[1], vertice7[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice4[0], vertice4[1], vertice4[2], vertice8[0], vertice8[1], vertice8[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice5[0], vertice5[1], vertice5[2], vertice6[0], vertice6[1], vertice6[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice5[0], vertice5[1], vertice5[2], vertice7[0], vertice7[1], vertice7[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice6[0], vertice6[1], vertice6[2], vertice8[0], vertice8[1], vertice8[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   CYLINDER, vertice7[0], vertice7[1], vertice7[2], vertice8[0], vertice8[1], vertice8[2], 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
   ]

cmd.load_cgo(obj,'box')