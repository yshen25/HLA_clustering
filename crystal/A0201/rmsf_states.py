#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell (rlc1@queensu.ca)

import math,sys
from pymol import cmd, stored

def distx2(x1,x2):
  """
  Calculate the square of the distance between two coordinates.
  Returns a float
  """
  distx2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
  return distx2


def rmsf_states(selection,byres=0,reference_state=1,debug=0):
  """
  AUTHOR

    Robert L. Campbell

  USAGE

    rmsf_states selection [,byres=0], [reference_state=1]

    Calculate the RMS fluctuation for each residue in a set of models
    (within a multi-state object).

    If you specify byres=1, then it will only calculate the RMSF for alpha-carbons (i.e.
    it modifies the selection by adding "and name CA"), but will modify the B-factors for
    all atoms in each residue to be the same as the alpha-carbons.

    If you use a selection that only includes only some atoms, only those atoms will have
    their B-factors modified (even with byres=1 specified).

    By default, the RMSF is calculated with respect to the first state, but you can
    change this by setting reference_state to the number of another state.
  """

# setting byres to true only makes sense if the selection is only for the alpha-carbons
  byres = int(byres)
  reference_state = int(reference_state)
  debug = int(debug)

# initialize the arrays used
  models = []
  coord_array = []
  chain=[]
  resi=[]
  resn=[]
  name=[]

  num_states = cmd.count_states(selection)

  if debug:
    print("There are %d states in your selection" % num_states)

# get models into models array
  for i in range(num_states):
    models.append(cmd.get_model(selection,state=i+1))

# extract coordinates and atom identification information out of the models
# loop over the states
  for i in range(num_states):
    coord_array.append([])
    chain.append([])
    resi.append([])
    resn.append([])
    name.append([])

    if byres:
# loop over the atoms in each state
      for j in range(len(models[0].atom)):
        if debug:
          print("i, j:", i,j)
        atom = models[i].atom[j]
        if atom.name == 'CA':
          coord_array[i].append(atom.coord)
          chain[i].append(atom.chain)
          resi[i].append(atom.resi)
          resn[i].append(atom.resn)
          name[i].append(atom.name)
    else:
# loop over the atoms in each state
      for j in range(len(models[0].atom)):
        if debug:
          print("i, j:", i,j)
        atom = models[i].atom[j]
        coord_array[i].append(atom.coord)
        chain[i].append(atom.chain)
        resi[i].append(atom.resi)
        resn[i].append(atom.resn)
        name[i].append(atom.name)

# initialize array
  diff2 = []

# calculate the square of the fluctuation
# = the sum of the distance between the atoms in each state from the reference state, squared
  for j in range(len(coord_array[0])):
    diff2.append(0.)
    for i in range(num_states):

# reference_state is 1-based, but coord_array is zero-based, so subtract 1 from 
# reference_state to get the correct reference_state
      diff2[j] += distx2(coord_array[i][j],coord_array[reference_state - 1][j])

# divide the fluctuation squared by the number of states and take the square root of that to get the RMSF
    diff2[j] = math.sqrt(diff2[j]/num_states)

# reset all the B-factors to zero
  cmd.alter(selection,"b=0")

  b_dict = {}

  sum_rmsf = 0
# if byres, alter all B-factors in a residue to be the same
  if byres:
    def b_lookup(chain, resi, name):
      if chain in b_dict:
        b = b_dict[chain][resi]
      else:
        b = b_dict[''][resi]
      return b
    stored.b = b_lookup

    for i in range(len(diff2)):
      sum_rmsf += diff2[i]
      b_dict.setdefault(chain[0][i], {})[resi[0][i]] = diff2[i]*8*math.pi**2
    mean_rmsf = sum_rmsf/len(diff2)
    cmd.alter(selection,'%s=stored.b(chain,resi,name)' % ('b'))

# if not byres, alter all B-factors individually according to the RMSF for each atom
  else:
    quiet=0
    def b_lookup(chain, resi, name):
      if chain in b_dict:
        b = b_dict[chain][resi][name]
      else:
        b = b_dict[''][resi][name]
      return b
    stored.b = b_lookup

    for i in range(len(diff2)):
      sum_rmsf += diff2[i]
      try:
        b_dict.setdefault(chain[0][i], {}).setdefault(resi[0][i], {})[name[0][i]] = diff2[i]*8*math.pi**2
      except KeyError:
        print(chain[0][i],resi[0][i],name[0][i])
    cmd.alter(selection,'%s=stored.b(chain,resi,name)' % ('b'))
    mean_rmsf = sum_rmsf/len(diff2)

  print("Mean RMSF for selection: %s = %g" % (selection, mean_rmsf))

cmd.extend("rmsf_states",rmsf_states)
