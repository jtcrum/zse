'''
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
'''

__all__ = ['get_orings','get_trings']

from zse.collections import get_ring_sizes
from zse.ring_utilities import *

# get_orings

def get_orings(atoms,index,code):
    '''
    Function to find all the rings asssociated with an oxygen atom in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')

    OUTPUTS:
    Class:      (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    # get possible rings, and max ring size
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)

    # get the closest neighbor of the oxygen, and find all possible rings
    # between that oxygen and its neighbor
    paths = get_paths(G,index,ring_sizes)

    # now we want to remove all the non ring paths
    paths = remove_non_rings(large_atoms, paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [int(len(p)/2) for p in paths]
    paths = [x for _,x in sorted(zip(ring_list,paths),reverse=True)]
    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms,paths)

    return ring_list, paths, ring_atoms

# get_trings
def get_trings(atoms,index,code):
    '''
    Function to find all the rings asssociated with a T-site in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')

    OUTPUTS:
    Class:      (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    #get possible rings, and max rings size
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)

    # to find all the rings associated with a T site, we need all the rings
    # associated with each oxygen bound to that T site. We will use networkx
    # neighbors to find those oxygens
    import networkx as nx
    paths = []
    for n in nx.neighbors(G,index):
        paths = paths+get_paths(G,n,ring_sizes)

    # now we want to remove all the non ring paths
    paths = remove_non_rings(large_atoms, paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [int(len(p)/2) for p in paths]
    paths = [x for _,x in sorted(zip(ring_list,paths),reverse=True)]
    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms,paths)

    return ring_list, paths, ring_atoms

# get_fwrings
