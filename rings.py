'''
The goal of this module is to identify the rings an oxygen is associated with if
a zeolite. This will give you some information as to the voids surrounding that
oxygen.

This module is in development form, and is not guarenteed to work
depending on the zeolite framework you use. It definitely works on
CHA, or potentially any zeolite with only one unique T Site.
'''

import networkx as nx
from ase.io import read, write
from ase import neighborlist
import numpy as np

def get_rings(atoms, index):

    '''
    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    '''

    cell = atoms.get_cell_lengths_and_angles()
    repeat = []

    for i,c in enumerate(cell):
        if c/2 < 8:
            repeat.append(2)
        else:
            repeat.append(1)
    atoms = atoms.repeat(repeat)
    center = atoms.get_center_of_mass()
    trans = center - atoms.positions[index]
    atoms.translate(trans)
    atoms.wrap()

    cutoff = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    for n in neighbs:
        fe = [n]
    fe.append(index)

    G.remove_edge(fe[0],fe[1])
    Class = []
    while len(Class)<6:
        try:
            path = nx.shortest_path(G,fe[0],fe[1])
        except:
            break
        Class.append(int(len(path)/2))
        for i in range(len(path)-3):
            G.remove_edge(path[i+1],path[i+2])
        Class.sort(reverse=True)
    return Class
