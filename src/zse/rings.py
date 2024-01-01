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

from zse.ring_utilities import (
    atoms_to_graph,
    paths_to_atoms,
    remove_dups,
    remove_geometric_dups,
    vertex_order,
)
from zse.ring_validation import crum, goetzke, sastre, vertex
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
