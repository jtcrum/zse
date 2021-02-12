__all__ = ['atoms_to_graph','get_paths','remove_non_rings','paths_to_atoms',
           'remove_dups','remove_labeled_dups','dict_to_atoms']

'''
This module contains utilities to be used by the rings.py module.
'''

import numpy as np
import math

def atoms_to_graph(atoms,index,max_ring):
    '''
    Helper function to repeat a unit cell enough times to capture the largest
    possible ring, and turn the new larger cell into a graph object.

    RETURNS:
    G = graph object representing zeolite framework in new larger cell
    large_atoms = ASE atoms object of the new larger cell framework
    repeat = array showing the number of times the cell was repeated: [x,y,z]
    '''

    # first, repeat cell, center the cell, and wrap the atoms back into the cell
    cell = atoms.get_cell_lengths_and_angles()[:3]
    repeat = []
    for i,c in enumerate(cell):
        if c/2 < max_ring/2+5:
            l = c
            re = 1
            while l/2 < max_ring/2+5:
                re += 1
                l = c*re

            repeat.append(re)
        else:
            repeat.append(1)
    large_atoms = atoms.copy()
    large_atoms = large_atoms.repeat(repeat)
    center = large_atoms.get_center_of_mass()
    trans = center - large_atoms.positions[index]
    large_atoms.translate(trans)
    large_atoms.wrap()

    # we need a neighborlist (connectivity matrix) before we can make our graph
    from ase import neighborlist
    cutoff = neighborlist.natural_cutoffs(large_atoms, mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(large_atoms)
    matrix = nl.get_connectivity_matrix(sparse=False)

    # now we make the graph
    import networkx as nx
    G = nx.from_numpy_matrix(matrix)

    return G, large_atoms, repeat

def get_paths(G,index,ring_sizes):
    '''
    Get all the paths (rings) between the index atom and its neighbor
    '''

    # first find the neighbor
    import networkx as nx
    neighbors = [n for n in nx.neighbors(G,index)]
    neighbor = neighbors[0]

    # next find all paths connecting index and neighbor
    paths = []
    for path in nx.all_simple_paths(G,index,neighbor,cutoff=max(ring_sizes)-1):
        if len(path) in ring_sizes:
            paths.append(path)

    return paths

def remove_non_rings(atoms, paths):
    # turn each path into a circular array to find and remove any duplicates
    paths = remove_dups(paths)

    # remove paths that contain smaller rings becuase these are likely not
    # actual rings
    # there are a ton of rules I have written to find these "secondary rings"
    # there is probably a better way to do this, I would love other's input.

    '''
    This portion doesn't seem to be necessary anymore.
    paths = remove_sec(paths)
    '''
    # this part is a trick to cut out some of the paths that are just
    # random walks through the framework, and not actual rings
    # here we have some cutoff distance for cross ring T-T distance
    # if the distance is less than the cutoff it probably isn't a ring
    # for odd MR we need a different function so there is no cross ring T-T


    delete = []
    for j,r in enumerate(paths):
        n = int(len(r)/2)
        cutoff = n/5
        if cutoff < 2:
            cutoff = 2
        if n%2 == 0 and n > 5:
            distances = []
            inner_flag = False
            for x in range(1,n,2):
                dist = atoms.get_distance(r[x],r[x+n])
                distances.append(dist)
                if dist < n-cutoff:
                    delete.append(j)
                    inner_flag = True
                    break
            if inner_flag == False:
                outer_flag = False
                for d in distances:
                    if d > n - math.floor(n/6):
                        outer_flag = True
                        break
                if outer_flag == False:
                    delete.append(j)
            if inner_flag == False and outer_flag == True:
                last_flag = False
                for d in distances:
                    if d > .94*n:
                        last_flag = False
                        break
                    else:
                        last_flag = True
                if last_flag:
                    delete.append(j)
        if n%2 != 0 and n > 5:
            r2 = r.copy()
            r2.append(r[:2])
            for x in range(1,n,2):
                dist1 = atoms.get_distance(r2[x],r2[x+n-1])
                dist2 = atoms.get_distance(r2[x],r2[x+n+1])
                if dist1 < n-cutoff or dist2 < n-cutoff:
                    delete.append(j)
                    break

    tmp_paths = paths.copy()
    paths = []
    for j,r in enumerate(tmp_paths):
        if j not in delete:
            paths.append(r)

    return paths

def remove_dups(paths):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    for i in range(len(paths)):
        for j in range((i+1), len(paths)):
            if i != j:
                st1 = set(paths[i])
                st2 = set(paths[j])
                if st1 == st2:
                    d.append(int(j))
    tmp_paths = []
    for i in range(len(paths)):
        if i not in d:
            tmp_paths.append(paths[i])
    paths = tmp_paths
    return paths

def remove_sec(paths):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    count2 = np.zeros(len(paths))

    for i in range(len(paths)):
        for j in range(i+1,len(paths)):
            if i!= j:
                ringi = paths[i]
                ringj = paths[j]
                ni = len(ringi)
                nj = len(ringj)
                if ni > nj and ni >= 16 and nj > 6:
                    count=0
                    for rj in ringj:
                        if rj in ringi:
                            count+=1
                    if count == nj/2:
                        count2[i]+=1
                    elif count > nj/2:
                        count2[i]+=2
                if nj >ni and nj >= 16 and ni > 6:
                    count=0
                    for ri in ringi:
                        if ri in ringj:
                            count+=1
                    if count == ni/2:
                        count2[j]+=1
                    elif count > ni/2:
                        count2[j]+=2
                if ni > nj and nj in [6,8]:
                    count=0
                    for rj in ringj:
                        if rj in ringi:
                            count+=1
                    if count >= nj-2:
                        count2[i]+=2
                if nj > ni and ni in [6,8]:
                    count=0
                    for ri in ringi:
                        if ri in ringj:
                            count+=1
                    if count >= ni-2:
                        count2[j]+=2
    for i,c in enumerate(count2):
        if c >=2:
            d.append(i)
    tmp_paths = []
    for i in range(len(paths)):
        if i not in d:
            tmp_paths.append(paths[i])
    paths = tmp_paths
    return paths

def paths_to_atoms(atoms,paths):
    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms if atom.index not in keepers]
    tmp_atoms = atoms.copy()
    del tmp_atoms[d]

    return tmp_atoms

def remove_labeled_dups(index_paths,label_paths,ring_sizes):

    # first make dictionaries
    label_rings = {}
    index_rings = {}
    for i,r in enumerate(label_paths):
        length = int(len(r)/2)
        if length not in label_rings:
            label_rings[length] = [r]
            index_rings[length] = [index_paths[i]]
        else:
            label_rings[length].append(r)
            index_rings[length].append(index_paths[i])

    # now we will remove duplicates
    for length in ring_sizes:
        ring_tlist = label_rings[length]
        ring_full  = index_rings[length]
        d = []
        for i in range(len(ring_tlist)):
            for j in range((i+1), len(ring_tlist)):
                st1 = ' '.join(map(str,ring_tlist[i]))
                st2 = ' '.join(map(str,ring_tlist[j]))
                st2_2 = ' '.join(map(str,reversed(ring_tlist[j])))
                if st2 in st1 + ' ' + st1 or st2_2 in st1 + ' ' + st1:
                    d.append(int(j))
        tmp1 = []
        tmp2 = []
        for i in range(len(ring_tlist)):
            if i not in d:
                tmp1.append(ring_tlist[i])
                tmp2.append(ring_full[i])
        label_rings[length] = tmp1
        index_rings[length] = tmp2

    return  index_rings, label_rings

def dict_to_atoms(index_paths,atoms):
    trajectories = {}
    com = atoms.get_center_of_mass()
    for length in index_paths.keys():
        paths = index_paths[length]
        tmp_traj = []
        for p in paths:
            tmp_atoms = paths_to_atoms(atoms,[p])
            trans = com-tmp_atoms[0].position
            tmp_atoms.translate(trans)
            tmp_atoms.wrap()
            tmp_traj+=[tmp_atoms]
        trajectories[length]=tmp_traj

    return trajectories
