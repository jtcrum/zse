'''
This is where I write new functions to test
before incorporating into the main code.
'''

__all__ = ['get_orings','get_trings','get_fwrings','get_geometry']

from zse.collections import get_ring_sizes, framework
from zse.ring_utilities import *
from zse.utilities import *
import numpy as np
from ase.geometry import get_distances
from zse.ring_validation import *

# get_orings

def get_orings(atoms,index,max_ring=12,validation='goetzke',cutoff=3.15):
    '''
    Function to find all the rings asssociated with an oxygen atom in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')
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
    # get possible rings, and max ring size
    ring_sizes = np.arange(3,max_ring+1)*2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)
    index = [atom.index for atom in large_atoms if atom.tag==index][0]


    # get the closest neighbor of the oxygen, and find all possible rings
    # between that oxygen and its neighbor
    if validation == 'goetzke':
        paths = goetzke(G,index,max_ring)
    else:
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

# get_trings
def get_trings(atoms,index,code,validation='cross_distance',cutoff=3.15):
    '''
    Function to find all the rings asssociated with a T-site in a zeolite
    framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (integer) index of the atom that you want to classify
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')
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


def get_geometry(atoms,index,code):
    from ase.geometry import get_distances
    from collections import defaultdict
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
    atoms = atoms.repeat(repeat)
    import networkx as nx
    com_dists = defaultdict(list)
    for n in nx.neighbors(G,index):
        c,paths,a = get_orings(large_atoms, n, code,validation = 'sastre')

        atoms,trans = center(atoms,large_atoms[n].tag)
        for p in paths:
            l = int(len(p)/2)
            ring_atoms = atoms[p]
            rp = ring_atoms.get_positions()
            com = rp.mean(axis=0)
            positions = atoms.get_positions()
            distances = get_distances(com,positions)[1][0]
            for i,d in enumerate(distances):
                if i not in p:
                    com_dists[l].append(d)
    return com_dists
