"""This module helps to identify the rings associated with an oxygen or
T-site in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, it's important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
"""

import warnings

from ase import Atoms

from zse.ring_utilities import (
    atoms_to_graph,
    paths_to_atoms,
    remove_dups,
    remove_geometric_dups,
    vertex_order,
)
from zse.ring_validation import crum, goetzke, sastre, vertex
from zse.utilities import center

__all__ = [
    "get_ordered_vertex",
    "get_rings",
    "get_unique_rings",
]


def get_rings(
    atoms: Atoms, index: int, validation: str | None = None, max_ring: int = 12
) -> tuple[list[int], list[list[int]], list[Atoms], Atoms]:
    """Find all the rings asssociated with an O-site or T-site in a
    zeolite framework.

    Args:
        atoms (Atoms): the zeolite framework to be analyzed, works best if you remove
            any adsorbates first
        index (int): index of the atom that you want to classify
        validation (str | None): Method with which to determin valid rings. Options are:
            None: No validation, return all rings found
            cross_distance: uses cross ring Si-Si distances
            d2: ensures each ring can't be decomposed into two smaller rings
            sphere: Checks that no non ring atoms are within some cutoff radius of
                    the center of mass of the ring.
                    Cutoff input is required for this method
            sp: Custom shortest path method, not very reliable
            vertex: Only valid for T-sites, ensures each ring contains two
                oxygen atoms bound to the T-site
            goetzke: Uses the algorithm presented by Goetzke and Klein to
                find all the rings up to a certain cutoff size for an atom.
                https://doi.org/10.1016/0022-3093(91)90145-V
        max_ring (int): Maximum size ring to search for. Time to compute scales with the
            maximum ring size.

    Returns:
        ring_list (list): The size of rings associated with the oxygen.
        paths (list): The actual atom indices that compose found rings.
        ring_atoms (list[Atoms]): all the rings found.
        atoms (Atoms): ASE Atoms object of the framework repeated large enough to
            visualize all the rings.
    """

    # need to multiply by two to consider oxygens
    max_ring *= 2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    graph, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = next(atom.index for atom in large_atoms if atom.tag == index)
    index_symbol = large_atoms[index].symbol

    # find cycles that don't contain any shortcuts
    paths = goetzke(graph, index, max_ring)

    # remove some cycles based on other validation rules
    if validation == "vertex":
        if index_symbol == "O":
            warnings.warn("Can't find vertex symbols of oxygen atoms", stacklevel=2)
            return False, False, False, False
        paths = vertex(paths)
    if validation == "sastre":
        paths = sastre(graph, paths, index_symbol)
    if validation == "crum":
        paths = crum(graph, paths, index_symbol)

    # convert the indices of the paths back to standard cell indices
    ring_list = [int(len(p) / 2) for p in paths]
    tmp_paths = [x for _, x in sorted(zip(ring_list, paths, strict=False))]
    paths = []
    for p in tmp_paths:
        temp = [large_atoms[i].tag for i in p]
        paths.append(temp)

    ring_list.sort()

    # make a collection of atoms objects to view the rings
    ring_atoms = [paths_to_atoms(large_atoms, [p]) for p in tmp_paths]
    atoms = atoms.repeat(repeat)

    return ring_list, paths, ring_atoms, atoms


def get_unique_rings(
    atoms: Atoms, tsites: list[int], validation: str | None = None, max_ring: int = 12
) -> tuple[list[int], list[list[int]], list[Atoms], Atoms]:
    """Find all the unique rings in a zeolite framework.
    This is accomplished by finding all the rings for each T-site, and then
    using T-O connectivity and geometry to remove duplicate rings.

    Args:
        atoms (Atoms): The ASE atoms object of the framework to classify
        tsites (list[int]): Indices of the unique tsites in the framework
            i.e. [101] for CHA (CHA has only one unique t-site)

    Returns:
        ring_list (list): List of the ring sizes found
        paths (list): List of the indices making each ring
        ring_atoms (list[Atoms]): Contains one image for each ring
        atoms (Atoms): ASE Atoms object of the framework repeated large enough
            to visualize all the rings.
    """

    paths = []
    for t in tsites:
        _c, r, _ra, a = get_rings(atoms, t, validation=validation, max_ring=max_ring)
        paths += r
    paths = remove_dups(paths)

    # remove duplicate rings based on geometry
    paths = remove_geometric_dups(a, paths)

    # sort rings from smallest to largest
    ring_list = [int(len(p) / 2) for p in paths]
    paths = [x for _, x in sorted(zip(ring_list, paths, strict=False))]
    ring_list.sort()

    ring_atoms = [paths_to_atoms(a, [p]) for p in paths]

    return ring_list, paths, ring_atoms, a


def get_ordered_vertex(
    atoms: Atoms, index: int, max_ring: int = 12
) -> tuple[str, list[list[int]], list[Atoms], Atoms]:
    """Find the vertex symbol of a given T-site, and return that
    vertex symbol with the rings listed in a specific order. This can be used
    to differientiate T-sites with the same vertex symbol, but different ring
    orientation around the T-site.

    Example:
        The vertex symbol of MOR T3 and MON T1 are both: 4·5_2·5·8_2·5·8_2
        This function would return the following however:
            MOR T3: 8_2•8_2•4•5_2•5•5
            MON T1: 8_2•8_2•5_2•4•5•5
        Letting us know that the connectivity of the rings around the T-site are
        different between these two T-sites. Visual inspection of the framework
        would tell us these two vertex symbols are different.

    Args:
        atoms (Atoms): The ASE atoms object of the framework to be analyzed.
            Works best if you remove any adsorbates first.
        index (int): Index of the atom that you want to classify.
            Must be a T-site and not an oxygen.
        max_ring (int): Maximum size ring to search for. Time to compute
            scales with the maximum ring size.

    Returns:
        ordered_vertex (str): The ordered vertex symbol for the T-site.
        paths (list): The actual atom indices that compose found rings.
            In the order they are presented in the ordered_vertex.
        ring_atoms (list[Atoms]): List of Atoms objects for all the rings found.
        atoms (Atoms): The ASE atoms object of the framework repeated large enough to
            visualize all the rings.
    """

    # need to multiply by two to consider oxygens
    max_ring *= 2

    # repeat the unit cell so it is large enough to capture the max ring size
    # also turn this new larger unit cell into a graph
    graph, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)
    index = next(atom.index for atom in large_atoms if atom.tag == index)
    index_symbol = large_atoms[index].symbol

    # find cycles that don't contain any shortcuts
    paths = goetzke(graph, index, max_ring)

    # remove some cycles based on other validation rules
    if index_symbol == "O":
        print("WARNING: Can't find vertex symbols of oxygen atoms")
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
    atoms, _vect = center(atoms, keep[0])
    atoms = atoms[keep]

    return ordered_vertex, paths, ring_atoms, atoms
