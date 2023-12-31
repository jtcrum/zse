"""
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
"""
from __future__ import annotations

import warnings

from zse.collections.framework import get_framework, get_ring_sizes
from zse.ring_utilities import (
    all_paths,
    atoms_to_graph,
    dict_to_atoms,
    get_paths,
    get_vertices,
    paths_to_atoms,
    remove_dups,
    remove_geometric_dups,
    remove_labeled_dups,
    shortest_valid_path,
    vertex_order,
)
from zse.ring_validation import (
    cross_distance,
    crum,
    goetzke,
    sastre,
    sp,
    sphere,
    vertex,
)
from zse.utilities import center, get_osites, site_labels


def get_rings(atoms, index, validation=None, max_ring=12):
    """
    Function to find all the rings asssociated with an O-site or T-site in a
    zeolite framework.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (int) index of the atom that you want to classify
    validation: (str) Method in which to determin valid rings.
                goetzke:    Uses the algorithm presented by Goetzke and Klein to
                            find all the rings up to a certain cutoff size for
                            an atom.
                            https://doi.org/10.1016/0022-3093(91)90145-V
    max_ring:   (int)       Maximum size ring to search for. Time to compute
                            scales with the maximum ring size.

    OUTPUTS:
    ring_list:  (list) The size of rings associated with the oxygen.
    paths:      (2d list) The actual atom indices that compose found rings.
    ring_atoms: (list of atoms objects) all the rings found.
    atoms:      (ase atoms object) of the framework repeated large enough to
                visualize all the rings.
    """

    # need to multiply by two to consider oxygens
    max_ring *= 2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = [atom.index for atom in large_atoms if atom.tag == index][0]
    index_symbol = large_atoms[index].symbol

    # find cycles that don't contain any shortcuts
    paths = goetzke(G, index, max_ring)

    if validation == "vertex":
        if index_symbol == "O":
            warnings.warn("WARNING: Can't find vertex symbols of oxygen atoms")
            return False, False, False, False
        paths = vertex(paths)
    if validation == "sastre":
        paths = sastre(G, paths, index_symbol)
    if validation == "crum":
        paths = crum(G, paths, index_symbol)

    # convert the indices of the paths back to standard cell indices
    ring_list = [len(p) // 2 for p in paths]
    tmp_paths = [x for _, x in sorted(zip(ring_list, paths))]
    paths = []
    for p in tmp_paths:
        temp = [large_atoms[i].tag for i in p]
        paths.append(temp)

    ring_list.sort()

    # make a collection of atoms objects to view the rings
    ring_atoms = [paths_to_atoms(large_atoms, [p]) for p in tmp_paths]
    atoms = atoms.repeat(repeat)

    return ring_list, paths, ring_atoms, atoms


def get_unique_rings(atoms, tsites, validation=None, max_ring=12):
    """
    Function to find all the unique rings in a zeolite framework.
    This is accomplished by finding all the rings for each T-site, and then
    using T-O connectivity and geometry to remove duplicate rings.

    INPUTS:
    atoms:          The ASE atoms object of the framework to classify
    tsites:         (list) Indices of the unique tsites in the framework
                    i.e. [101] for CHA (CHA has only one unique t-site)

    OUTPUTS:
    ring_list:      (list) List of the ring sizes found
    paths:          (list) list of the indices making each ring
    ring_atoms:     (list of ASE atoms) Contains one image for each ring
    atoms:          (ASE atoms object) of the framework repeated large enough to
                    visualize all the rings.
    """

    paths = []
    for t in tsites:
        _, r, _, a = get_rings(atoms, t, validation=validation, max_ring=max_ring)
        paths += r
    paths = remove_dups(paths)

    # remove duplicate rings based on geometry
    paths = remove_geometric_dups(a, paths)

    # sort rings from smallest to largest
    ring_list = [len(p) // 2 for p in paths]
    paths = [x for _, x in sorted(zip(ring_list, paths))]
    ring_list.sort()

    ring_atoms = [paths_to_atoms(a, [p]) for p in paths]

    return ring_list, paths, ring_atoms, a


def get_ordered_vertex(atoms, index, max_ring=12):
    """
    Function to find the vertex symbol of a given T-site, and return that
    vertex symbol with the rings listed in a specific order. This can be used
    to differientiate T-sites with the same vertex symbol, but different ring
    orientation around the T-site.

    i.e.
        The vertex symbol of MOR T3 and MON T1 are both: 4·5_2·5·8_2·5·8_2
        This function would return the following however:
            MOR T3: 8_2•8_2•4•5_2•5•5
            MON T1: 8_2•8_2•5_2•4•5•5
        Letting us know that the connectivity of the rings around the T-site are
        different between these two T-sites. Visual inspection of the framework
        would tell us these two vertex symbols are different.

    INPUTS:
    atoms:      (ASE atoms object) the zeolite framework to be analyzed
                works best if you remove any adsorbates first
    index:      (int) Index of the atom that you want to classify.
                      Must but a T-site and not an oxygen.
    max_ring:   (int) Maximum size ring to search for. Time to compute
                      scales with the maximum ring size.

    OUTPUTS:
    ordered_vertex:  (str) The ordered vertex symbol for the T-site.
    paths:      (2d list) The actual atom indices that compose found rings.
                In the order they are presented in the ordered_vertex.
    ring_atoms: (list of atoms objects) for all the rings found.
    atoms:      (ase atoms object) of the framework repeated large enough to
                visualize all the rings.
    """

    # need to multiply by two to consider oxygens
    max_ring *= 2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = [atom.index for atom in large_atoms if atom.tag == index][0]
    index_symbol = large_atoms[index].symbol

    # find cycles that don't contain any shortcuts
    paths = goetzke(G, index, max_ring)

    # remove some cycles based on other validation rules
    if index_symbol == "O":
        warnings.warn("WARNING: Can't find vertex symbols of oxygen atoms")
        return False, False, False, False
    paths = vertex(paths)

    # convert the indices of the paths back to standard cell indices
    ordered_vertex, paths = vertex_order(paths)
    tmp_paths = paths
    paths = []
    for p in tmp_paths:
        temp = [large_atoms[i].tag for i in p]
        paths.append(temp)

    # get the ordered vertex symbol

    # make a collection of atoms objects to view the rings
    ring_atoms = [paths_to_atoms(large_atoms, [p]) for p in tmp_paths]
    atoms = atoms.repeat(repeat)
    keep = []
    for path in paths:
        for x in path:
            if x not in keep:
                keep.append(x)
    atoms, _ = center(atoms, keep[0])
    atoms = atoms[keep]

    return ordered_vertex, paths, ring_atoms, atoms


""" DEPRECRATED FUNCTIONS """


# get_orings
def get_orings(atoms, index, code, validation="cross_distance", cutoff=3.15):
    """
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
    """
    # get possible rings, and max ring size
    ring_sizes = get_ring_sizes(code) * 2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, _ = atoms_to_graph(atoms, index, max_ring)
    index = [atom.index for atom in large_atoms if atom.tag == index][0]

    # get the closest neighbor of the oxygen, and find all possible rings
    # between that oxygen and its neighbor
    paths = get_paths(G, index, ring_sizes)
    paths = remove_dups(paths)

    # now we want to remove all the non ring paths
    # the validation method will determine which set of rules to use
    # to eliminate non ring paths
    if validation == "sp":
        paths = sp(G, paths)
    elif validation == "d2":
        raise ValueError("d2 validation not implemented")
    elif validation == "sphere":
        if cutoff is None:
            raise ValueError(
                "INPUT ERROR: Validation with geometry requires cutoff in Å, however, cutoff not set."
            )
        paths = sphere(large_atoms, paths, cutoff)
    elif validation == "cross_distance":
        paths = cross_distance(large_atoms, paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [len(p) // 2 for p in paths]
    tmp_paths = [x for _, x in sorted(zip(ring_list, paths), reverse=True)]
    paths = []
    for p in tmp_paths:
        temp = [large_atoms[i].tag for i in p]
        paths.append(temp)

    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms, tmp_paths)

    return ring_list, paths, ring_atoms


# get_trings
def get_trings(atoms, index, code, validation="cross_distance", cutoff=3.15):
    """
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
    """
    # get possible rings, and max rings size
    ring_sizes = get_ring_sizes(code) * 2
    max_ring = max(ring_sizes)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = [atom.index for atom in large_atoms if atom.tag == index][0]

    # to find all the rings associated with a T site, we need all the rings
    # associated with each oxygen bound to that T site. We will use networkx
    # neighbors to find those oxygens
    import networkx as nx

    paths = []
    for n in nx.neighbors(G, index):
        paths = paths + get_paths(G, n, ring_sizes)

    # Since we found the rings for each oxygen attached to the T-site,
    # there will be duplicate rings. Let's remove those.
    paths = remove_dups(paths)

    # now we want to remove all the non ring paths, the method for determining
    # valid rings is designated with the validation input variable.
    if validation == "sp":
        paths = sp(G, paths)
    elif validation == "d2":
        raise ValueError("d2 validation not implemented")
    elif validation == "sphere":
        if cutoff is None:
            raise ValueError(
                "INPUT ERROR: Validation with geometry requires cutoff in Å, however, cutoff not set."
            )
        paths = sphere(large_atoms, paths, cutoff)
    elif validation == "cross_distance":
        paths = cross_distance(large_atoms, paths)

    # finally organize all outputs: list of ring sizes, atom indices that make
    # ring paths, and an atoms object that shows all those rings
    ring_list = [len(p) // 2 for p in paths]
    paths2 = [x for _, x in sorted(zip(ring_list, paths), reverse=True)]
    tmp_paths = [x for _, x in sorted(zip(ring_list, paths), reverse=True)]
    paths = []
    for p in tmp_paths:
        temp = [large_atoms[i].tag for i in p]
        paths.append(temp)
    ring_list.sort(reverse=True)

    ring_atoms = paths_to_atoms(large_atoms, paths2)

    return ring_list, paths, ring_atoms


# get_fwrings
def get_fwrings(code, validation="cross_distance", cutoff=3.15):
    """
    Function to find all the unique rings in a zeolite framework.

    INPUTS:
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')

    OUTPUTS:
    index_paths:    (dictionary) {ring length : indices of atoms in ring}
    label_paths:    (dictionary) {ring length : site labels of atoms in ring}
    trajectories:   (dictionary of atoms objects)  -
                    {ring length : [atoms objects for that size ring]}
    """

    # First, get some basic info we will need about the framework
    # atoms object, possible ring sizes, and the index of each o-site type
    atoms = get_framework(code)
    ring_sizes = get_ring_sizes(code)
    _, _, oinds = get_osites(code)

    # now we find all the rings associated with each oxygen type in the fw
    # this in theory should find us every possible type of ring in the fw
    paths = []
    for o in oinds:
        paths += get_orings(atoms, o, code, validation=validation, cutoff=cutoff)[1]
    paths = remove_dups(paths)

    # now we want to get the t-site and o-site labels
    repeat = atoms_to_graph(atoms, 0, max(ring_sizes) * 2)[2]
    atoms = atoms.repeat(repeat)
    labels = site_labels(atoms, code)

    # now convert all the paths from index form into site label form
    index_paths = paths
    label_paths = []
    for path in index_paths:
        l = [labels[p] for p in path]
        label_paths.append(l)

    # now we want to remove duplicate rings based on the label_paths
    # this will also give us the paths in a conveinent dictionary
    # we will check the labels and geometry to remove duplicates
    index_paths, label_paths = remove_labeled_dups(
        index_paths, label_paths, ring_sizes, atoms
    )

    # last but not least, let's make an atoms object for each ring type so that
    # we can visualize them later
    # this will be a dictionary of trajectories
    trajectories = dict_to_atoms(index_paths, atoms)

    return index_paths, label_paths, trajectories


# get_vertex_symbols
def get_vertex_symbols(code, index):
    """
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
    """

    # get the possible rings of this framework from the IZA
    # the maximum ring size determines how many times to repeat the unit cell
    ring_sizes = get_ring_sizes(code) * 2
    max_ring = max(ring_sizes)

    # get an atoms object of the framework
    atoms = get_framework(code)

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    G, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = [atom.index for atom in large_atoms if atom.tag == index][0]

    # get the T-site and O-site labels for each atom in the framework
    atoms = atoms.repeat(repeat)
    labels = site_labels(atoms, code)

    # get each o-o pair for the T-site
    vertices = get_vertices(G, index)

    # set up a dictionary to stare results
    vertex_symbols = {}

    # got through each o-o pair and find the shortest paths
    traj = []
    for i, v in enumerate(vertices):
        o1 = v[0]
        o2 = v[1]
        v_label = f"{labels[large_atoms[o1].tag]}-{labels[large_atoms[o2].tag]}"

        # this finds the shortest path between the two
        _, l = shortest_valid_path(G, o1, o2, index)
        # this finds all valid paths of that length between the two
        paths = list(all_paths(G, o1, o2, index, l))
        traj += [paths_to_atoms(large_atoms, paths)]
        tmp_paths = []
        for p in paths:
            temp = [large_atoms[x].tag for x in p]
            tmp_paths.append(temp)
        vertex_symbols[f"{i + 1}:{v_label}"] = list(tmp_paths)

    return vertex_symbols, traj
