# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:23:13 2020

@author: nicme
"""

# 2D structure --> LAMMPS input file

node_list = "node_list.dat"
bond_list = "link_info.dat"

import numpy as np
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
Setup
'''

def nodes_and_bonds(node_filename, bond_filename,
                    scale_down_factor, frac):

    #import bonds
    bonds = []
    with open(bond_filename, 'r') as node_data:
        lines = node_data.readlines()
        for line in lines:
            line_list = [int(x) for x in line.split(',')]      
            bonds.append(line_list)
    node_data.close()
    #return bonds
    
    #import nodes and scale down
    nodes = []
    with open(node_filename, 'r') as node_data:
        lines = node_data.readlines()
        for line in lines:
            line_list = [float(x)/scale_down_factor 
                         for x in line.split(',')]           
            nodes.append((line_list))
    
    node_data.close()
    
    # make bonds list with coords instead of node numbers
    
    bonds_coord = []
    for bond in bonds:
        node_1 = int(bond[0]) - 1
        node_2 = int(bond[1]) - 1
        
        coord_1 = nodes[node_1]
        coord_2 = nodes[node_2]
    
        bond_coord = [coord_1, coord_2]
        bonds_coord.append(bond_coord)
    
    #cut nodes, then bonds
    # try to make node coordinates smaller

    nodes_array = np.array(nodes)
    xmax = max(nodes_array[:,0])*frac
    ymax = max(nodes_array[:,1])*frac
    zmax = max(nodes_array[:,2])*frac
    # convert back to nested list   
    
    short_nodes = []
    for node in nodes_array:   
        node = list(node)
        if node[0] <= xmax and node[1] <= ymax and node[2] <= zmax:
            short_nodes.append(node)        
            
    # go through bonds_coord 
    # and remove missing nodes from short nodes list
    short_bonds_coord = []
    for bond in bonds_coord:
        if bond[0] in short_nodes and bond[1] in short_nodes:
                short_bonds_coord.append(bond)
    
    # revert from bonds_coord to bond index
    short_bonds = []
    for bond_coord in short_bonds_coord:
        node_1_coord = bond_coord[0] 
        node_2_coord = bond_coord[1] 
        
        node_1 = short_nodes.index(node_1_coord) + 1
        node_2 = short_nodes.index(node_2_coord) + 1
        
        bond = [node_1, node_2]
        short_bonds.append(bond)
        
    return short_nodes, short_bonds

def plot_mesh(nodes, bonds):

    #nodes = np.array(nodes)
    
    X_vals = np.array(nodes)[:, 0]
    Y_vals = np.array(nodes)[:, 1]
    Z_vals = np.array(nodes)[:, 2]

    ax = plt.axes(projection='3d')
    ax.scatter3D(X_vals, Y_vals, Z_vals, color = 'green', s=2.5)
    
    
    for bond in np.array(bonds):
        node_1 = bond[0] - 1
        node_2 = bond[1] - 1
        
        x1, y1, z1 = nodes[node_1]
        x2, y2, z2 = nodes[node_2]        
        ax.plot([x1, x2], [y1, y2], [z1, z2], color = 'maroon', linewidth = 3)
        
#    plt.title('ECM Mesh')
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.zlabel('z')
    plt.show()

'''
ECM mesh gemoetry calculations
'''

def connection_lengths(nodes, bonds):
    connection_lengths = []
    
    for i in range(len(bonds)):
        connection = bonds[i]
        
        node_1 = connection[0]-1
        node_2 = connection[1]-1 # because python lists begin with a zero
        
        X_node_1 = nodes[node_1][0]
        Y_node_1 = nodes[node_1][1]
        Z_node_1 = nodes[node_1][2]
        
        X_node_2 = nodes[node_2][0]
        Y_node_2 = nodes[node_2][1]
        Z_node_2 = nodes[node_2][2]

        length = np.sqrt((X_node_2 - X_node_1)**2 + 
                         (Y_node_2 - Y_node_1)**2 +
                         (Z_node_2 - Z_node_1)**2)

        connection_lengths.append(length)
    
    return connection_lengths

def angles(nodes, bonds):
    node_angle_triplet = [] # triplet of connected nodes
    node_angle = []
    
    
    bond_combinations = list(itertools.combinations(bonds, 2))
    #bond_combinations_section = bond_combinations[0:5] # for testing
    #print(bond_combinations_section)
    # look at senior design notes for visualization #
    
    
    for combo in bond_combinations:
        bond_1 = combo[0] #[node_1, node_2]
        bond_2 = combo[1] #[node_3, node_4]
        
        # find node in common
        # common node = value = c
        c = [value for value in bond_1 if value in bond_2]
        #print(c)

        if c == []:
            pass
        else:
            c = c[0]
            triplet = bond_1 + bond_2
            # remove c and add to the end
            triplet.remove(c)
            triplet.remove(c) 
            triplet.append(c)

            node_angle_triplet.append(triplet)        
            
    node_angle_triplet = list(node_angle_triplet)
    # node_angle_triplet works
        
    for triplet in node_angle_triplet:
        node_1_label = triplet[0] - 1 # that's how Python labels
        node_2_label = triplet[1] - 1
        node_3_label = triplet[2] - 1 # middle point
        
        node_1 = nodes[node_1_label]
        node_2 = nodes[node_2_label]
        node_3 = nodes[node_3_label]
        
        d_13 = np.asarray(node_1) - np.asarray(node_3)
        d_23 = np.asarray(node_2) - np.asarray(node_3)
        
        angle_cos = np.dot(d_13, d_23)/(np.linalg.norm(d_13) * np.linalg.norm(d_23))
        angle = np.arccos(angle_cos)
        angle = np.degrees(angle)
        
        
        node_angle.append(angle)

    return node_angle_triplet, node_angle
    
def uniform_mass_cmds(uniform_mass):
    masses_1 = "1   "+str(uniform_mass)
    masses_2 = "2   "+str(uniform_mass)
    masses = [masses_1, masses_2]
    
    return masses

def harm_bond_cmds(harm_bond_coeff, bond_lengths):
    bond_cmds = []
    
    for i in range(len(bond_lengths)):
        bond_cmd = str(i+1) + '  '
        bond_cmd = bond_cmd + str(harm_bond_coeff) + '  '
        bond_cmd = bond_cmd + str(bond_lengths[i])
        bond_cmds.append(bond_cmd)
        
    return bond_cmds

def angle_coeffs_cmds(angle_coeff, node_angle):
    angle_cmds = [] #angle commands
    for i in range(len(node_angle)):
        triplet_angle = node_angle[i]
        
        angle_cmd = str(i+1)+'  '+str(angle_coeff)
        angle_cmd = angle_cmd + '  ' + str(triplet_angle)
        
        angle_cmds.append(angle_cmd)
        
    return angle_cmds
    
def atom_cmds(nodes):
    atom_list = []
    for i in range(len(nodes)): 
        # first 3 terms
        string = str(i+1) +'   0  ' + str(1) # --> Masses type 1
        # x and y coordinates
        string = string +'  ' + str(nodes[i][0]) + '  ' + str(nodes[i][1])
        # z coordinate
        string = string +'  ' + str(nodes[i][2]) # change to str(nodes[i][2]) in 3D
        atom_list.append(string)
        
    return atom_list

def bond_cmds(bonds):
    bond_cmds = []
    for i in range(len(bonds)):
        bond = bonds[i]
        node_1 = bond[0]
        node_2 = bond[1]
        bond_cmd = str(i+1) + '  ' + str(i+1) # --> from Bond_coeffs
        bond_cmd = bond_cmd + '  ' + str(node_1)
        bond_cmd = bond_cmd + '  ' + str(node_2)
        
        bond_cmds.append(bond_cmd)
        
    return bond_cmds
    
def angle_cmds(node_triplet):
    angle_cmds = []
    
    for i in range(len(node_triplet)):
        
        triplet = node_triplet[i]
        node_1 = triplet[0]
        node_2 = triplet[2] # connecting node is last in triplet format
        node_3 = triplet[1]
        
        angle_cmd = str(i+1) + '  ' + str(i+1) # same order as Bond Coeffs
        angle_cmd = angle_cmd + '  ' + str(node_1)
        angle_cmd = angle_cmd + '  ' + str(node_2)
        angle_cmd = angle_cmd + '  ' + str(node_3)
        
        angle_cmds.append(angle_cmd)
    
    return angle_cmds

def write_input_file(nodes, bonds, node_triplet,
                     Masses, Bond_Coeffs, Angle_Coeffs,
                     Atoms, Bonds, Angles, 
                     filename):

    open(filename, 'w').close() # erases previous content of in_file
    in_file = open(filename, 'w')
    
    
    # first line
    L1 = ['LAMMPS 3D ECM input file\n\n']
    in_file.writelines(L1)
    
    # atoms, bonds, angles
    L2 = [str(len(nodes)) + '  atoms\n',
          str(len(bonds)) + '  bonds\n',
          str(len(node_triplet)) + '  angles\n\n']
    in_file.writelines(L2)
    
    # atom types, bond types, angle types
    L3 = [str(len(Masses)) + '  atom types\n',
          str(len(Bond_Coeffs)) + '  bond types\n',
          str(len(Angle_Coeffs)), '  angle types\n\n']
    in_file.writelines(L3)
    
    # xlo xhi ylo yhi
    nodes_array = np.array(nodes)
    xlo = min(nodes_array[:,0]) # use zero instead
    xhi = 1.1 * max(nodes_array[:,0])
    ylo = min(nodes_array[:,1]) # use zero instead
    yhi = 1.1 * max(nodes_array[:,1])
    zlo = min(nodes_array[:,2]) # use zero instead
    zhi = 1.1 * max(nodes_array[:,2])
    
    L4 = [str(0) + ' ' + str(xhi) + ' xlo' + ' xhi\n',
          str(0) + ' ' + str(yhi) + ' ylo' + ' yhi\n',
          str(0) + ' ' + str(zhi) + ' zlo' + ' zhi\n\n']

    in_file.writelines(L4)
    
    
    # Masses
    L5 = ['Masses\n\n']
    for mass in Masses:
        line = mass + '\n'
        L5.append(line)
    L5.append('\n')
    in_file.writelines(L5)
    
    # Bond Coeffs
    L6 = ['Bond Coeffs\n\n']
    for bond_coeff in Bond_Coeffs:
        line = bond_coeff + '\n'
        L6.append(line)
    L6.append('\n')
    in_file.writelines(L6)
    
    # Angle Coeffs
    L7 = ['Angle Coeffs\n\n']
    for angle_coeff in Angle_Coeffs:
        line = angle_coeff + '\n'
        L7.append(line)
    L7.append('\n')
    in_file.writelines(L7)
    
    # Atoms
    L8 = ['Atoms\n\n']
    for atom in Atoms:
        line = atom + '\n'
        L8.append(line)
    L8.append('\n')
    in_file.writelines(L8)
    
    # Bonds
    L9 = ['Bonds\n\n']
    for bond in Bonds:
        line = bond + '\n'
        L9.append(line)
    L9.append('\n')
    in_file.writelines(L9)
    
    # Angles
    L10 = ['Angles\n\n']
    for angle in Angles:
        line = angle + '\n'
        L10.append(line)
    L10.append('\n')
    in_file.writelines(L10)
    
    in_file.close()
