'''
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
'''

import networkx as nx
from ase.io import read, write
from ase import neighborlist
import numpy as np
import math
from zse import substitute

def get_rings(atoms, index):

    '''
    WARNING: This is old and does not work for all framework types.
    WARNING: Use the updated get_orings function below instead.


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

    cutoff = neighborlist.natural_cutoffs(atoms,mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    for n in neighbs:
        if atoms[n].symbol == 'Si':
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

def get_orings(atoms, index,possible):

    '''
    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    possible: (list) of the types of rings known to be present in the zeolite
              framework you are studying. This information is available on IZA
              or in the collections module of this package.
    Returns: Class - The size of rings associated with the oxygen.
             Rings - The actual atom indices that compose those rings.
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

    cutoff = neighborlist.natural_cutoffs(atoms,mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    for n in neighbs:
        if atoms[n].symbol == 'Si':
            fe = [n]
    fe.append(index)

    G.remove_edge(fe[0],fe[1])
    tmpClass = []
    rings = []
    while len(tmpClass)<6:
        try:
            path = nx.shortest_path(G,fe[0],fe[1])
        except:
            break
        if len(path) in possible:
            tmpClass.append(int(len(path)/2))
            rings.append(path)
        length = len(path)
        for n in np.arange(math.ceil(length/4),math.ceil(length - length/4)):
            G.remove_edge(path[n],path[n+1])

    rings = remove_dups(rings)
    rings = remove_sec(rings)
    Class = []
    for r in rings:
        Class.append(int(len(r)/2))
    Class.sort(reverse=True)
    return Class, rings

def find_o_rings(G,index,possible):
    '''
    This is a helper function for the get_trings function.
    It won't do much on its own.
    '''
    rings = []
    neighbs = nx.neighbors(G,index)
    oxygen = []
    for n in neighbs:
        oxygen.append(n)
    for i in range(4):
        G2 = G.copy()
        oneighbs = nx.neighbors(G2,oxygen[i])
        neighbors = []
        for o in oneighbs:
            neighbors.append(o)
        neighbor = neighbors[0]
        G2.remove_edge(neighbor, oxygen[i])
        tmp_class = []

        while len(tmp_class)<6:
            try:
                path = nx.shortest_path(G2,oxygen[i],neighbor)
            except:
                break
            if len(path) in possible:
                tmp_class.append(int(len(path)/2))
                rings.append(path)
            length = len(path)
            for n in np.arange(math.ceil(length/4),math.ceil(length - length/4)):
                G2.remove_edge(path[n],path[n+1])

    return rings

def remove_dups(rings):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    for i in range(len(rings)):
        for j in range((i+1), len(rings)):
            if i != j:
                st1 = set(rings[i])
                st2 = set(rings[j])
                if st1 == st2:
                    d.append(int(j))
    paths = []
    for i in range(len(rings)):
        if i not in d:
            paths.append(rings[i])
    return paths

def remove_sec(rings):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            if i!= j:
                ringi = rings[i]
                ringj = rings[j]
                ni = len(ringi)
                nj = len(ringj)
                if ni > nj:
                    count=0
                    for rj in ringj:
                        if rj in ringi:
                            count+=1
                    if count > nj/2:
                        d.append(i)
                if nj >ni:
                    count=0
                    for ri in ringi:
                        if ri in ringj:
                            count+=1
                    if count > ni/2:
                        d.append(j)

    paths = []
    for i in range(len(rings)):
        if i not in d:
            paths.append(rings[i])
    return paths

def get_trings(atoms,index,possible):
    '''
    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    possible: (list) of the types of rings known to be present in the zeolite
              framework you are studying. This information is available on IZA
              or in the collections module of this package.
    Returns: Rings - The actual atom indices that compose those rings.
             atoms2 - An atoms object with the desired T Site changed to an
             Aluminum atom (just for visual purposes), and all atoms removed
             except for those that share a ring with the T Site provided.
    '''

    atoms2 = atoms.copy()
    cell = atoms2.get_cell_lengths_and_angles()
    repeat = []

    for i,c in enumerate(cell):
        if c/2 <12:
            l = c
            re = 2
            while l/2 < 12:
                l = c*re
                re+=1
            repeat.append(re)
        else:
            repeat.append(1)
    atoms2 = atoms2.repeat(repeat)
    center = atoms2.get_center_of_mass()
    trans = center - atoms2.positions[index]
    atoms2.translate(trans)
    atoms2.wrap()

    cutoff = neighborlist.natural_cutoffs(atoms2, mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms2)
    matrix = nl.get_connectivity_matrix(sparse = False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)
    rings = find_o_rings(G,index,possible)
    paths = remove_dups(rings)
    paths = remove_sec(paths)
#     print('Unique Rings')
#     for r in paths:
#         print(int(len(r)/2),r)
    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms2 if atom.index not in keepers]
    atoms2 = substitute.tsub(atoms2,index,'Al')
    del atoms2[d]

    return paths, atoms2
