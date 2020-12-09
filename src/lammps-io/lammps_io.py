'''
Filename: lammps_io.py
Programmed by: Niccolo Meniconi and Mike Bernard
Date: 2020-12-08
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from time import time
from tqdm import tqdm


def make_graph_from_file(node_filepath, bond_filepath, scale_factor=1.0, timeit=False):
    ''' Read in node and bond data from CSV files.
    :param node_filename: `str` the file path to a CSV file of node_1,node_2 data
    :param bond_filename: `str` the file path to a CSV file of x,y or x,y,z data for each node
    :param scale_factor: `float` a factor by which to scale the coordinates of each node
    :return: `networkx.Graph` a graph of the node and edge data
    '''
    # TODO: implement BOUND_FACTOR to ignore nodes which fall outside the bounding
    # box of all the nodes scaled by some factor

    print('Reading data in... ', end='')

    if timeit:
        start = time()
    
    G = nx.Graph()

    # read in bonds from file
    # bond number = line number - 1 due to zero-indexing
    # list of [node_1, node_2] representing the numbers of the two nodes a bond connects
    bonds = np.array([[int(num) - 1 for num in line.strip().split(',')]
                        for line in open(bond_filepath, 'r')])
    
    # read in nodes from file
    # list of [x, y] or [x, y, z] representing coordinates of each node
    nodes = np.array([[float(coord)*scale_factor for coord in line.strip().split(',')]
                      for line in open(node_filepath, 'r')])

    # add all nodes to the graph
    G.add_nodes_from([
        (i, {'coords': np.array(nodes[i])})
        for i in range(np.shape(nodes)[0])
    ])

    # add all edges to the graph
    for (node_i, node_j) in bonds:
        G.add_edge(node_i, node_j)

    if timeit:
        print(time() - start)
    
    print('Done')
    return G


def calc_connection_lengths(G, timeit=False):
    ''' Calculate the lengths of all the bonds in the graph.
    :param G: `networkx.Graph` the graph of nodes and edges
    :return: `dict` {node_i: [node_j: dist_ij, ..., node_m: dist_im]}
    '''
    print('Calculating connection lengths... ', end='')

    if timeit:
        start = time()

    how_edgy = {}
    
    for node in G.nodes:
        how_edgy[node] = {}
        coords_1 = G.nodes[node]['coords']
        conns = G.adj[node]
        for flight_of_the_conn in conns:
            coords_2 = G.nodes[flight_of_the_conn]['coords']  # heh
            how_edgy[node][flight_of_the_conn] = np.linalg.norm(coords_2 - coords_1)
    
    if timeit:
        print(time() - start)

    print('Done')
    return how_edgy


def find_triplets(G, timeit=False):
    ''' Find all triplets of nodes connected by edges in a graph.
    :param G: `networkx.Graph` the graph of nodes and edges
    :return: `np.ndarray` of [node_1, node_2, node_3], sorted by node
             number, containing only unique values.
    '''
    print('Finding triplets... ', end='')

    if timeit:
        start = time()
    
    triplets = []
    for node_0 in G.nodes:
        for node_1 in G.adj[node_0]:
            for node_2 in G.adj[node_1]:
                if len(set([node_0, node_1, node_2])) == 3:
                    triplets.append([node_0, node_1, node_2])
    
    triplets = np.array(triplets)
    triplets.sort()
    
    if timeit:
        print(time() - start)
    
    print('Done')
    return triplets


def calc_angles(G, triplets, timeit=False):
    ''' Calculate the angles formed by all triplets of nodes in a graph.
    :param G: `networkx.graph` the graph of nodes and edges
    :param triplets: `(n, 3) numpy.ndarray` array of triplets of nodes
                     connected by edges, sorted by node number, containing
                     only unique triplets, as output by find_triplets(G)
    :return: `dict` {node_i: {node_j: {node_k: angle_ijk, ..., node_z: angle_ijz}}}
             where angle_ijn is the angle formed by the connection between node_i
             and node_j and the connection between node_j and node_n; the smaller
             of the two possible angles is given, in radians, as `float`
    '''
    print('Calculating angles... ')
    if timeit:
        start = time()
    
    angles = {}

    for (i, j, k) in tqdm(triplets):
        if i not in angles.keys():
            angles[i] = {}
        if j not in angles[i].keys():
            angles[i][j] = {}
        
        v_ji = G.nodes[i]['coords'] - G.nodes[j]['coords']
        v_jk = G.nodes[k]['coords'] - G.nodes[j]['coords']

        angle = np.arccos(np.dot(v_ji, v_jk)/(np.linalg.norm(v_ji)*np.linalg.norm(v_jk)))
        angles[i][j][k] = angle

    if timeit:
        print(time() - start)

    print('Done')
    return angles


# TODO: implement write_lammps_input_file
# def write_lammps_input_file(G, writepath, mass, harmonic_bond_coeff, angle_coeff,
#                             connection_lengths, angles, node_coords, bonds, triplets,
#                             timeit=False):
#     ''' Writes an input file for LAMMPS. '''
#     print(f'Writing LAMMPS input file to {writepath}... ', end='')

#     if timeit:
#         start = time()

#     mass_cmds = [f'1   {mass}', f'2   {mass}']
#     harmonic_bonds_cmds = [
#         f'{i}  {harmonic_bond_coeff}  {connection_lengths[i]}'
#                       for i in range(len(connection_lengths))
#     ]
#     angle_coeffs_cmds = [
#         f'{i+1}  {angle_coeff}  {angles[i]}'
#         for i in angles)
#     ]
#     atom_cmds = [
#         f'{i+1}   1  {node_coords[i][0]}  {str(node_coords[i])[1:-1]}'
#         for i in range(len(node_coords))
#     ]
#     bond_cmds = [
#         f'{i+1}  {i+1}  {bonds[i][0]}  {bonds[i][1]}'
#         for i in range(len(bonds))
#     ]
#     angle_cmds = [
#         f'{i+1}  {i+1}  {i}  {j}  {k}'
#         for (i, j, k) in triplets
#     ]

#     if timeit:
#         print(time() - start)
    
#     print('Done')


def show_graph(G):
    ''' Display the graph using Matplotlib.
    :param G: `networkx.Graph` the graph of nodes and edges
    :return: None
    '''
    start = time()
    node_coords = np.array([coords for coords in nx.get_node_attributes(G, 'coords').values()])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot nodes
    ax.scatter(node_coords[:,0], node_coords[:,1], node_coords[:,2])
    
    # plot edges
    edge_coords = []
    for node in G.nodes:
        coords_1 = G.nodes[node]['coords']
        conns = G.adj[node]
        for conn in conns:
            coords_2 = G.nodes[conn]['coords']
            edge_coords.append([coords_1, coords_2])
    edge_coords = np.array(edge_coords)
    
    for node_1, node_2 in edge_coords:
        ax.plot([node_1[0], node_2[0]], [node_1[1], node_2[1]], zs=[node_1[2], node_2[2]])

    print(time() - start)
    plt.show()


def main():
    BOND_STYLE = 'harmonic'  # TODO: add description
    ANGLE_STYLE = 'harmonic'  # TODO: add description

    UNIFORM_MASS = 1e-11  # TODO: add units, description
    HARMONIC_BOND_COEFFICIENT = 1e-11  # TODO: add units, description
    ANGLE_COEFFICIENT = 0.0  # TODO: add units, description

    SCALE_FACTOR = 1.0  # scale coordinates of each node by this much
    BOUND_FACTOR = 0.8  # ignore nodes outside the bounding box scaled by this much

    nodefile_path = "/Users/mike/Downloads/MD_Dec2020/ECM_3D/node_list.dat"
    bondfile_path = "/Users/mike/Downloads/MD_Dec2020/ECM_3D/link_info.dat"
    outfile_path = '/Users/mike/Downloads/MD_Dec2020/ECM_3D/3D_input_file.data'

    G = make_graph_from_file(nodefile_path, bondfile_path)
    node_coords = [G.nodes[i]['coords'] for i in G.nodes]
    edges = [e for e in G.edges]
    connection_lengths = calc_connection_lengths(G)
    triplets = find_triplets(G)
    angles = calc_angles(G, triplets)
    # show_graph(G)

    # write_lammps_input_file(
    #     G, outfile_path, UNIFORM_MASS, HARMONIC_BOND_COEFFICIENT, ANGLE_COEFFICIENT,
    #     connection_lengths, angles, node_coords, edges, triplets
    # )


if __name__ == '__main__':
    main()
