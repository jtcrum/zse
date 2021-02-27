'''
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
'''

__all__ = ['get_orings','get_trings','get_fwrings']

from zse.collections import get_ring_sizes, framework
from zse.ring_utilities import *
from zse.utilities import *
import numpy as np
from ase.geometry import get_distances

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
    ring_list:      (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    # get possible rings, and max ring size
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph2(atoms,index,max_ring)
    index = [atom.index for atom in large_atoms if atom.tag==index][0]


    # get the closest neighbor of the oxygen, and find all possible rings
    # between that oxygen and its neighbor
    paths = get_paths(G,index,ring_sizes)

    # now we want to remove all the non ring paths
    paths = remove_non_rings(large_atoms, paths)

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
    ring_list:      (list) The size of rings associated with the oxygen.
    paths:      (2d array) The actual atom indices that compose found rings.
    ring_atoms: (ASE atoms object) all the rings found
    '''
    #get possible rings, and max rings size
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

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

    # now we want to remove all the non ring paths
    #paths = remove_non_rings(large_atoms, paths)

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

# get_fwrings
def get_fwrings(code):
    '''
    Function to find all the unique rings in a zeolite framework.

    INPUTS:
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')

    OUTPUTS:
    index_paths:    (dictionary) {ring length : indices of atoms in ring}
    label_paths:    (dictionary) {ring length : site labels of atoms in ring}
    trajectories:   (dictionary of atoms objects)  -
                    {ring length : [atoms objects for that size ring]}
    '''

    # First, get some basic info we will need about the framework
    # atoms object, possible ring sizes, and the index of each o-site type
    atoms = framework(code)
    ring_sizes = get_ring_sizes(code)
    osites,omults,oinds = get_osites(code)

    # now we find all the rings associated with each oxygen type in the fw
    # this in theory should find us every possible type of ring in the fw
    paths = []
    for o in oinds:
        paths += get_orings(atoms,o,code)[1]
    paths = remove_dups(paths)

    # now we want to get the t-site and o-site labels
    repeat = atoms_to_graph(atoms,0,max(ring_sizes)*2)[2]
    atoms = atoms.repeat(repeat)
    labels = site_labels(atoms,code)

    # now convert all the paths from index form into site label form
    index_paths = paths
    label_paths = []
    for path in index_paths:
        l = []
        for p in path:
            l.append(labels[p])
        label_paths.append(l)

    # now we want to remove duplicate rings based on the label_paths
    # this will also give us the paths in a conveinent dictionary
    index_paths, label_paths = remove_labeled_dups(index_paths, label_paths, ring_sizes)

    # last but not least, let's make an atoms object for each ring type so that
    # we can visualize them later
    # this will be a dictionary of trajectories
    trajectories = dict_to_atoms(index_paths,atoms)

    return index_paths, label_paths, trajectories

def atoms_to_graph2(atoms,index,max_ring,scale = True):
    '''
    Helper function to repeat a unit cell enough times to capture the largest
    possible ring, and turn the new larger cell into a graph object.

    RETURNS:
    G = graph object representing zeolite framework in new larger cell
    large_atoms = ASE atoms object of the new larger cell framework
    repeat = array showing the number of times the cell was repeated: [x,y,z]
    '''

    # first scale the unit cell so the average Si-Si distance = 3.1 Å
    if scale:
        atoms = scale_cell(atoms)


    # repeat cell, center the cell, and wrap the atoms back into the cell
    cell = atoms.cell.cellpar()[:3]
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

    # remove atoms that won't contribute to wrings
    from ase.geometry import get_distances
    cell = large_atoms.get_cell()
    pbc = [1,1,1]
    p1 = large_atoms[index].position
    positions = large_atoms.get_positions()
    distances = get_distances(p1,positions)[1][0]

    delete = []
    for i,l in enumerate(distances):
        if l>max_ring/2+5:
            delete.append(i)
    inds = [atom.index for atom in large_atoms]
    large_atoms.set_tags(inds)
    atoms = large_atoms.copy()
    del large_atoms[delete]

    matrix = np.zeros([len(large_atoms),len(large_atoms)]).astype(int)
    positions = large_atoms.get_positions()



    tsites = [atom.index for atom in large_atoms if atom.symbol != 'O']
    tpositions = positions[tsites]
    osites = [atom.index for atom in large_atoms if atom.index not in tsites]
    opositions = positions[osites]
    distances = get_distances(tpositions,opositions)[1]

    for i,t in enumerate(tsites):
        dists = distances[i]
        idx = np.nonzero(dists<2)[0]
        for o in idx:
            matrix[t,osites[o]]=1
            matrix[osites[o],t]=1
    # now we make the graph
    import networkx as nx
    G = nx.from_numpy_matrix(matrix)
    # G.remove_nodes_from(delete)
    return G, large_atoms, repeat
