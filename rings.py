'''
This module is in development form, and is not guarenteed to work
depending on the zeolite framework you use. It definitely works on
CHA, or potentially any zeolite with only one unique T Site. 
'''

import numpy as np
from ase import neighborlist
from ase.io import read, write

def zeoc(Adjacency, Source, LengthOfCycle, CurrentPath=[0], UsedNodes=[0],  iteration=0, CurrentVertex=0):
    global CycleList
    if iteration == 0:
        CycleList= []
        CurrentVertex = Source
        UsedNodes = Source
        CurrentPath = [Source]
        iteration += 1
        zeoc(Adjacency, Source, LengthOfCycle, CurrentPath, [UsedNodes], iteration, CurrentVertex)
    else:
        if iteration == 1:
            UsedNodes[0] = []
        iteration += 1
        if len(CurrentPath) == LengthOfCycle and CurrentPath[-1] == Source:
            CycleList.append(CurrentPath.copy())
            Cycle = CurrentPath.copy()
            Adjacency = zerout(CurrentPath, Adjacency)
        elif len(CurrentPath) < LengthOfCycle:
            CurrentNeighbors = [x for x in range(len(Adjacency[CurrentVertex,:])) if Adjacency[CurrentVertex,x] != 0 ]
            for i in range(len(CurrentNeighbors)):
                Used = False
                CurrentVertex = CurrentNeighbors[i]
                for j in range(len(UsedNodes)):
                    if CurrentNeighbors[i] == UsedNodes[j]:
                        Used = True
                if not Used:
                    CurrentPath.append(CurrentVertex)
                    UsedNodes.append(CurrentVertex)
                    NodeIndex = len(UsedNodes)-1
                    Cycles, Adjacency = zeoc(Adjacency, Source, LengthOfCycle, CurrentPath, UsedNodes, iteration, CurrentVertex)
                    UsedNodes[NodeIndex] = []
                    del CurrentPath[-1]
                    #CurrentPath[-1] = []
            mmm = 0
        elif len(CurrentPath) == LengthOfCycle and CurrentPath[-1] != Source:
            UsedNodes[-1] = []
            CurrentPath[-1] = []
    Cycles = CycleList
    return Cycles, Adjacency

def zerout(CurrentPath, Adjacency):
    NewAdjacency = Adjacency
    temp = CurrentPath.copy()
    del temp[0]
    del temp[0]
    del temp[-1]
    del temp[-1]
    GraphSize = max(Adjacency.shape)
    Z = np.zeros([1, GraphSize])
    Z1 = np.zeros([GraphSize,1])
    for i in range(len(temp)):
        NewAdjacency[temp[i],:] = 0
        NewAdjacency[:,temp[i]] = 0
    return NewAdjacency

def get_rings(atoms, index, possible_rings):
    # first we need the connectivity matrix
    # atoms: ase Atoms object of a zeolite structure
    # index: (integer) index of the oxygen atom that you want to classify
    # possible_rings: (array) what size rings are possible in this zeolite framework?
    # Check IZA if you don't know. Example: CHA has 8-,6-,4-MR i.e. [8,6,4]

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

    Class = []
    pr = possible_rings
    pr.sort()
    for i in pr:
        ring,m = zeoc(m,index,i*2+1)
        if ring:
            for r in ring:
                temp=r
                Class.append(int((len(temp)-1)/2))
    Class.sort(reverse=True)
    return Class
