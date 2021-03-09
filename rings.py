'''
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
'''

__all__ = ['get_orings','get_trings','get_fwrings','get_vertex_symbols']

from zse.collections import get_ring_sizes, framework
from zse.ring_utilities import *
from zse.utilities import *
from zse.ring_validation import *

# get_orings

def get_orings(atoms,index,code,validation='cross_distance',cutoff=3.15):
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
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

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
def get_fwrings(code,validation='cross_distance',cutoff=3.15):
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
        paths += get_orings(atoms,o,code, validation = validation, cutoff = cutoff)[1]
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
    # we will check the labels and geometry to remove duplicates
    index_paths, label_paths = remove_labeled_dups(index_paths, label_paths, ring_sizes,atoms)

    # last but not least, let's make an atoms object for each ring type so that
    # we can visualize them later
    # this will be a dictionary of trajectories
    trajectories = dict_to_atoms(index_paths,atoms)

    return index_paths, label_paths, trajectories

def get_vertex_symbols(code,index):
    '''
    Function to find all the the shortest rings connecting each
    oxygen-oxygen pair associated with a T-site

    INPUTS:
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')
    index:      (integer) index of the T-site that you want to classify

    OUTPUTS:
    vertex_symbols: (dictionary) keys are the O-site labels, and values are the
                    indices of the atoms connecting those oxygens
    trajectory:     (ASE atoms objects) that include the rings for each
                    oxygen-oxygen pair
    '''

    # get the possible rings of this framework from the IZA
    # the maximum ring size determines how many times to repeat the unit cell
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)

    # get an atoms object of the framework
    atoms = framework(code)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)
    index = [atom.index for atom in large_atoms if atom.tag==index][0]

    # get the T-site and O-site labels for each atom in the framework
    atoms = atoms.repeat(repeat)
    labels = site_labels(atoms,code)

    # get each o-o pair for the T-site
    vertices = get_vertices(G,index)

    # set up a dictionary to stare results
    vertex_symbols = {}

    # got through each o-o pair and find the shortest paths
    traj = []
    for i,v in enumerate(vertices):
        o1 = v[0]
        o2 = v[1]
        v_label = '{0}-{1}'.format(labels[large_atoms[o1].tag],labels[large_atoms[o2].tag])

        # this finds the shortest path between the two
        path, l = shortest_valid_path(G,o1,o2,index)
        # this finds all valid paths of that length between the two
        paths = [p for p in all_paths(G,o1,o2,index,l)]
        traj+=[paths_to_atoms(large_atoms,paths)]
        tmp_paths = []
        for p in paths:
            temp = []
            for x in p:
                temp.append(large_atoms[x].tag)
            tmp_paths.append(temp)
        vertex_symbols['{0}:{1}'.format(i+1,v_label)] = [path for path in tmp_paths]

    return vertex_symbols, traj
