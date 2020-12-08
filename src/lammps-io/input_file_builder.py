# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:53:25 2020

@author: nicme
"""

#from input_file_reader import * # gives me an annoing error
                                # import functions one by one
from input_file_reader import nodes_and_bonds
from input_file_reader import connection_lengths, angles
from input_file_reader import plot_mesh
from input_file_reader import uniform_mass_cmds
from input_file_reader import harm_bond_cmds
from input_file_reader import angle_coeffs_cmds
from input_file_reader import atom_cmds, bond_cmds
from input_file_reader import angle_cmds
from input_file_reader import write_input_file


#import itertools
import numpy as np
from time import time
start = time()

node_list = "/Users/mike/Downloads/MD_Dec2020/ECM_3D/node_list.dat"
bond_list = "/Users/mike/Downloads/MD_Dec2020/ECM_3D/link_info.dat"

'''
ECM parameter definition: uniform parameters
'''

bond_style  = 'harmonic'
angle_style = 'harmonic'

uniform_mass = 1e-11
harm_bond_coeff = 1e-11
#harm_bond_coeff = 0.0
angle_coeff = 0.0


filename = '/Users/mike/Downloads/MD_Dec2020/ECM_3D/3D_input_file.data' # input file name


'''
Scale down
'''

scale_factor = 1 # scale down coordinate points
frac = 0.8

## work on changing import_bond_list later

'''
Setup
'''

nodes, bonds = nodes_and_bonds(node_list, bond_list, 
                               scale_factor, frac)

## bugs
## graph has unconnected nodes


'''
### find way to download plotly

#plot_mesh(short_nodes, short_bonds)
'''


'''
#ECM mesh geometry calculations
'''

bond_lengths = connection_lengths(nodes, bonds)
    # bond_lengths[i] = bond length of bonds[i]

node_triplet, node_angle = angles(nodes, bonds)
    # node_angle[i] = angle between nodes 1->2->3 of node_triplet[i]
    # middle node: 3rd node

####################################################
# ECM parameter definition: non-uniform parameters #
######### ignore, work on in the future ############
####################################################


#mass_range = [1, 2, 0.1]
#bond_coeff_range = [1, 2, 0.1]
#angle_coeff_range = [1, 2, 0.1]


'''
#functions that writes data file strings 
''' 



Masses = uniform_mass_cmds(uniform_mass)
# string that go under Masses
# len(Masses) = atom types

Bond_Coeffs = harm_bond_cmds(harm_bond_coeff, bond_lengths)
# string that go under Bond Coeffs
# bond_cmds[n] corresponds to the comand inputting the bond[n]
# len(Bond_Coeffs) = bond types

Angle_Coeffs = angle_coeffs_cmds(angle_coeff, node_angle)
# strings that go under Angle Coeffs
# Angle_Coeffs[n] corresponds to the comand inputting the node_triplet[n]
# len(Angle_Coeffs) = angle types

Atoms = atom_cmds(nodes)
# strings that go under Atoms

Bonds = bond_cmds(bonds)
# strings that go under Bonds

Angles = angle_cmds(node_triplet)
# strings that go under Angles

'''
#Writing input file
'''

write_input_file(nodes, bonds, node_triplet,
                     Masses, Bond_Coeffs, Angle_Coeffs,
                     Atoms, Bonds, Angles, 
                     filename)
end = time()
dur = end - start
print(f'Total Time: {dur:.3f} sec')
