'''
This module will find all rings for a particular O-site, T-site, or framework
up to a cutoff set by the user. I recommend setting the cutoff to 12 or less
because computation time increases significantly with cutoff.
'''
__all__ = ['get_orings','get_trings']

from zse.ring_utilities import *
from zse.utilities import *
from zse.ring_validation import *

import numpy as np

def get_orings(atoms,index,max_ring = 12,validation='d2',cutoff=3.15):
    '''
    Function to find all the rings asssociated with an oxygen atom in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    max_ring:   (integer) the maximum number of T-sites allowed in a ring
    validation: (str) Method in which to determin valid rings.
                cross_distance: uses cross ring Si-Si distances
                d2:             ensures each ring can't be decomposed into two
                                smaller rings
                sphere:         Checks that no non ring atoms are within some
                                cutoff radius of the center of mass of the
                                ring.
                                Cutoff input is required for this method
                sp:             Custum shortest path method, not very reliable
    cutoff:     (float) Value required for the sphere validation method

    OUTPUTS:
    ring_list:  (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    ring_sizes = np.arange(3,max_ring+1)*2
    max_ring *=2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)
    index = [atom.index for atom in large_atoms if atom.tag==index][0]


    # get the closest neighbor of the oxygen, and find all possible rings
    # between that oxygen and its neighbor
    paths = get_paths(G,index,ring_sizes)
    paths = remove_dups(paths)

    # now we want to remove all the non ring paths
    # the validation method will determine which set of rules to use
    # to eliminate non ring paths
    if validation == 'sp':
        paths = sp(G,paths)
    if validation == 'd2':
        paths = d2(G,paths)
    if validation =='sphere':
        if cutoff == None:
            print('INPUT ERROR: Validation with geometry requires cutoff in Å, however, cutoff not set.')
            return
        paths = sphere(large_atoms,paths,cutoff)
    if validation == 'cross_distance':
        paths = cross_distance(large_atoms,paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [int(len(p)/2) for p in paths]
    tmp_paths = [x for _,x in sorted(zip(ring_list,paths),reverse=True)]
    paths = []
    for p in tmp_paths:
        temp = []
        for i in p:
            temp.append(large_atoms[i].tag)
        paths.append(temp)

    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms,tmp_paths)

    return ring_list, paths, ring_atoms

def get_trings(atoms,index,max_ring = 12,validation='d2',cutoff=3.15):
    '''
    Function to find all the rings asssociated with a T-site in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    max_ring:   (integer) the maximum number of T-sites allowed in a ring
    validation: (str) Method in which to determin valid rings.
                cross_distance: uses cross ring Si-Si distances
                d2:             ensures each ring can't be decomposed into two
                                smaller rings
                sphere:         Checks that no non ring atoms are within some
                                cutoff radius of the center of mass of the
                                ring.
                                Cutoff input is required for this method
                sp:             Custum shortest path method, not very reliable
    cutoff:     (float) Value required for the sphere validation method

    OUTPUTS:
    ring_list:      (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    #get possible rings, and max rings size
    ring_sizes = np.arange(3,max_ring+1)*2
    max_ring *=2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)
    index = [atom.index for atom in large_atoms if atom.tag==index][0]

    # to find all the rings associated with a T site, we need all the rings
    # associated with each oxygen bound to that T site. We will use networkx
    # neighbors to find those oxygens
    import networkx as nx
    paths = []
    for n in nx.neighbors(G,index):
        paths = paths+get_paths(G,n,ring_sizes)

    # Since we found the rings for each oxygen attached to the T-site,
    # there will be duplicate rings. Let's remove those.
    paths = remove_dups(paths)

    # now we want to remove all the non ring paths, the method for determining
    # valid rings is designated with the validation input variable.
    if validation == 'sp':
        paths = sp(G,paths)
    if validation == 'd2':
        paths = d2(G,paths)
    if validation =='sphere':
        if cutoff == None:
            print('INPUT ERROR: Validation with geometry requires cutoff in Å, however, cutoff not set.')
            return
        paths = sphere(large_atoms,paths,cutoff)
    if validation == 'cross_distance':
        paths = cross_distance(large_atoms,paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [int(len(p)/2) for p in paths]
    paths2 = [x for _,x in sorted(zip(ring_list,paths),reverse=True)]
    tmp_paths = [x for _,x in sorted(zip(ring_list,paths),reverse=True)]
    paths = []
    for p in tmp_paths:
        temp = []
        for i in p:
            temp.append(large_atoms[i].tag)
        paths.append(temp)
    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms,paths2)

    return ring_list, paths, ring_atoms
