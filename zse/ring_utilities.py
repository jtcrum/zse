"""This module contains utilities to be used by the rings.py module."""

from collections import defaultdict as dd
from itertools import permutations as perm

import networkx as nx
import numpy as np
from ase import Atoms
from ase.geometry import get_distances

from zse.utilities import center

__all__ = [
    "all_paths",
    "atoms_to_graph",
    "dict_to_atoms",
    "get_paths",
    "get_vertices",
    "is_valid",
    "paths_to_atoms",
    "remove_dups",
    "remove_geometric_dups",
    "remove_labeled_dups",
    "shortest_valid_path",
    "vertex_order",
]


def atoms_to_graph(atoms: Atoms, index: int, max_ring: int) -> tuple[nx.Graph, Atoms, list[int]]:
    """Repeat a unit cell enough times to capture the largest
    possible ring, and turn the new larger cell into a graph object.

    Args:
        atoms (Atoms): ASE atoms object of the zeolite framework
        index (int): index of the T atom to be analyzed
        max_ring (int): size of the largest ring to be analyzed

    Returns:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        large_atoms (Atoms): ASE atoms object of the new larger cell framework
        repeat (list[int]): array showing the number of times the cell was repeated: [x,y,z]
    """

    # repeat cell, center the cell, and wrap the atoms back into the cell
    cell = atoms.cell.cellpar()[:3]
    repeat = []
    for c in cell:
        if c / 2 < max_ring / 2 + 5:
            length = c
            re = 1
            while length / 2 < max_ring / 2 + 5:
                re += 1
                length = c * re

            repeat.append(re)
        else:
            repeat.append(1)
    large_atoms = atoms.copy()
    large_atoms = large_atoms.repeat(repeat)
    center = large_atoms.get_center_of_mass()
    trans = center - large_atoms.positions[index]
    large_atoms.translate(trans)
    large_atoms.wrap()

    cell = large_atoms.get_cell()
    p1 = large_atoms[index].position
    positions = large_atoms.get_positions()
    distances = get_distances(p1, positions)[1][0]

    delete = []
    for i, length in enumerate(distances):
        if length > max_ring / 2 + 5:
            delete.append(i)
    inds = [atom.index for atom in large_atoms]
    large_atoms.set_tags(inds)
    atoms = large_atoms.copy()
    del large_atoms[delete]

    matrix = np.zeros([len(large_atoms), len(large_atoms)]).astype(int)
    positions = large_atoms.get_positions()

    tsites = [atom.index for atom in large_atoms if atom.symbol != "O"]
    tpositions = positions[tsites]
    osites = [atom.index for atom in large_atoms if atom.index not in tsites]
    opositions = positions[osites]
    distances = get_distances(tpositions, opositions)[1]

    for i, t in enumerate(tsites):
        dists = distances[i]
        idx = np.nonzero(dists < 2)[0]
        for o in idx:
            matrix[t, osites[o]] = 1
            matrix[osites[o], t] = 1

    graph = nx.from_numpy_array(matrix)
    return graph, large_atoms, repeat


def remove_dups(paths: list) -> list:
    """A helper function for get_orings and get_trings to remove duplicate paths.

    Args:
        paths (list): list of lists containing the paths

    Returns:
        paths (list): list of lists containing the paths with duplicates removed
    """
    d = []
    for i in range(len(paths)):
        for j in range((i + 1), len(paths)):
            if i != j:
                st1 = set(paths[i])
                st2 = set(paths[j])
                if st1 == st2:
                    d.append(int(j))
    tmp_paths = [paths[i] for i in range(len(paths)) if i not in d]
    paths = tmp_paths
    return paths


def paths_to_atoms(atoms: Atoms, paths: list) -> Atoms:
    """Convert a list of paths (lists of indices) to an ASE atoms object.

    Args:
        atoms (Atoms): ASE atoms object of the zeolite framework
        paths (list): list of lists containing the paths

    Returns:
        tmp_atoms (Atoms): ASE atoms object containing only the atoms in the paths
    """
    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms if atom.index not in keepers]
    tmp_atoms = atoms.copy()
    del tmp_atoms[d]

    return tmp_atoms


def remove_geometric_dups(atoms: Atoms, paths: list) -> list:
    """A helper function for get_orings and get_trings to remove geometrically
    duplicate paths.

    Args:
        atoms (Atoms): ASE atoms object of the zeolite framework
        paths (list): list of lists containing the paths

    Returns:
        unique_paths (list): list of lists containing the paths with geometric duplicates removed
    """
    unique_paths = []
    unique_distances = []
    for p in paths:
        atoms, _trans = center(atoms, p[0])
        dists = []
        for x in range(len(p) - 1):
            for r in range(x + 1, len(p)):
                dist = round(atoms.get_distance(p[x], p[r]), 2)
                dists.append(dist)
        dists.sort()

        if dists not in unique_distances:
            unique_paths.append(p)
            unique_distances.append(dists)
    return unique_paths


def get_vertices(graph: nx.Graph, index: int) -> list:
    """Get all the vertices (pairs of neighbors) connected to the index atom.

    Args:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        index (int): index of the T atom to be analyzed

    Returns:
        vertices (list): list of lists containing the vertices (pairs of neighbors)
    """
    vertices = []
    neighbors = list(nx.neighbors(graph, index))
    length = len(neighbors)
    for j in range(length - 1):
        for k in range(j + 1, length):
            vertices.append([neighbors[j], neighbors[k]])  # noqa: PERF401
    return vertices


def shortest_valid_path(
    graph: nx.Graph, o1: int, o2: int, index: int, length: int
) -> tuple[list, int]:
    """Find the shortest valid path between two oxygen atoms that goes through
    the index atom.

    Args:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        o1 (int): index of the first oxygen atom
        o2 (int): index of the second oxygen atom
        index (int): index of the T atom to be analyzed
        l (int): maximum length of the path

    Returns:
        p (list): list containing the shortest valid path
        length (int): length of the shortest valid path
    """
    graph_ = graph.copy()
    graph_.remove_node(index)
    flag = True
    l2 = 6
    p_idx = []
    while flag:
        paths = nx.all_simple_paths(graph_, o1, o2, l2 - 1)
        for path_ in paths:
            p_idx.append(index)
            if len(path_) == l2:
                flag, _j = is_valid(graph, path_)
            if not flag:
                break
        if flag:
            l2 += 2
        if l2 > length:
            p_idx = [1]
            break
    return p_idx, len(p_idx)


def is_valid(graph: nx.Graph, path: list) -> tuple[bool, int]:
    """Check if a path is valid (i.e., does not have shortcuts).

    Args:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        path (list): list containing the path to be checked

    Returns:
        flag (bool): True if the path is not valid, False if it is valid
        j (int): index of the first node in the path that creates a shortcut
    """
    path_length = len(path)
    flag = False
    for j in range(1, path_length - 1, 2):
        for k in range(j + 2, path_length, 2):
            node = path[j]
            sp = nx.shortest_path(graph, node, path[k])

            if len(sp) < k - j + 1 and len(sp) < path_length - (k - j) + 1:
                flag = True
                break
        if flag:
            break
    return flag, j


def all_paths(graph: nx.Graph, o1: int, o2: int, index: int, path_length: int) -> list:
    """Find all valid paths between two oxygen atoms that go through
    the index atom.

    Args:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        o1 (int): index of the first oxygen atom
        o2 (int): index of the second oxygen atom
        index (int): index of the T atom to be analyzed
        path_length (int): length of the path

    Returns:
        all_paths (list): list of lists containing all valid paths
    """
    all_paths = []
    graph_ = graph.copy()
    graph_.remove_node(index)
    paths = nx.all_simple_paths(graph_, o1, o2, path_length)
    for path in paths:
        path.append(index)
        if len(path) == path_length:
            flag, _j = is_valid(graph, path)
            if not flag and path not in all_paths:
                all_paths.append(path)
    return all_paths


def vertex_order(all_paths: list) -> tuple[str, list]:
    """Order the paths based on the oxygen pairs and their weights.

    Args:
        all_paths (list): list of lists containing the paths

    Returns:
        ordered_v (str): string representing the ordered vertices
        new_r (list): list of lists containing the ordered paths
    """
    oxygens = []
    o_pair_paths = dd(list)
    o_pair_sizes = dd(int)
    o_pair_counts = dd(lambda: 0)
    o_pair_weights = dd(lambda: 0)

    for path in all_paths:
        if path[1] not in oxygens:
            oxygens.append(path[1])
        if path[-1] not in oxygens:
            oxygens.append(path[-1])
        oxys = sorted([str(path[1]), str(path[-1])])
        o_pair_paths["-".join(oxys)].append(path)
        o_pair_sizes["-".join(oxys)] = len(path)
        o_pair_counts["-".join(oxys)] += 1
        o_pair_weights["-".join(oxys)] += len(path)

    perms = list(perm(oxygens))
    weights = []
    order = [[0, 2], [0, 1], [1, 2], [2, 3], [3, 0], [1, 3]]
    for p in perms:
        w = []
        for o in order:
            try:
                k = "-".join(sorted([str(p[o[0]]), str(p[o[1]])]))
                w.append(o_pair_weights[k])
            except:  # noqa: E722, PERF203
                pass
        weights.append(w)

    zipped_lists = zip(weights, perms, strict=False)
    sp = sorted(zipped_lists, reverse=True)
    tuples = zip(*sp, strict=False)
    weights, perms = [list(tuple) for tuple in tuples]

    new_r = []
    counts = []
    sizes = []
    oxygens = perms[0]
    weights = weights[0]

    for o in order:
        try:
            k = "-".join(sorted([str(oxygens[o[0]]), str(oxygens[o[1]])]))
            sizes.append(o_pair_sizes[k])
            counts.append(o_pair_counts[k])
            for x in o_pair_paths[k]:
                new_r.append(x)  # noqa: PERF402
        except:  # noqa: E722, PERF203
            pass

    ordered_v = []
    for c, s in zip(counts, sizes, strict=False):
        if c == 1:
            ordered_v.append(f"{int(s / 2)}")
        elif c == 0:
            ordered_v.append("*")
        else:
            ordered_v.append(f"{int(s / 2)}_{c}")
    ordered_v = "•".join(ordered_v)
    c = [int(len(x) / 2) for x in new_r]

    return ordered_v, new_r


def dict_to_atoms(index_paths: dict, atoms: Atoms) -> dict[int, list[Atoms]]:
    """Convert a dictionary of paths (lists of indices) to a dictionary of ASE atoms objects.

    Args:
        index_paths (dict): dictionary where keys are path lengths and values are lists of paths
        atoms (Atoms): ASE atoms object of the zeolite framework

    Returns:
        trajectories (dict): dictionary where keys are path lengths and values are lists of
            ASE atoms objects
    """
    trajectories = {}
    com = atoms.get_center_of_mass()
    for length, paths in index_paths.items():
        tmp_traj = []
        for p in paths:
            tmp_atoms = paths_to_atoms(atoms, [p])
            trans = com - tmp_atoms[0].position
            tmp_atoms.translate(trans)
            tmp_atoms.wrap()
            tmp_traj += [tmp_atoms]
        trajectories[length] = tmp_traj

    return trajectories


def remove_labeled_dups(
    index_paths: dict, label_paths: list, ring_sizes: list, atoms: Atoms
) -> tuple[dict, dict]:
    """A helper function for get_orings and get_trings to remove duplicate paths based
    on their labels.

    Args:
        index_paths (dict): dictionary where keys are path lengths and values are lists of paths
        label_paths (list): list of lists containing the labels of the paths
        ring_sizes (list): list of integers representing the sizes of rings to be analyzed
        atoms (Atoms): ASE atoms object of the zeolite framework

    Returns:
        index_rings (dict): dictionary where keys are path lengths and values are lists of paths
        label_rings (dict): dictionary where keys are path lengths and values are lists of labels
    """
    # first make dictionaries
    label_rings = {}
    index_rings = {}
    for i, r in enumerate(label_paths):
        length = int(len(r) / 2)
        if length not in label_rings:
            label_rings[length] = [r]
            index_rings[length] = [index_paths[i]]
        else:
            label_rings[length].append(r)
            index_rings[length].append(index_paths[i])

    # now we will remove duplicates
    for length in ring_sizes:
        ring_tlist = label_rings[length]
        ring_full = index_rings[length]
        d = []
        for i in range(len(ring_tlist)):
            for j in range((i + 1), len(ring_tlist)):
                st1 = " ".join(map(str, ring_tlist[i]))
                st2 = " ".join(map(str, ring_tlist[j]))
                st2_2 = " ".join(map(str, reversed(ring_tlist[j])))
                if st2 in st1 + " " + st1 or st2_2 in st1 + " " + st1:
                    p = ring_full[i]
                    atoms, _trans = center(atoms, p[0])
                    cross1 = []
                    for x in range(1, len(p) + 1, 2):
                        for r in range(x + 2, len(p) + 1, 2):
                            dist = round(atoms.get_distance(p[x], p[r]), 1)
                            cross1.append(dist)
                    cross1.sort()

                    p = ring_full[j]
                    atoms, _trans = center(atoms, p[0])
                    cross2 = []
                    for x in range(1, len(p) + 1, 2):
                        for r in range(x + 2, len(p) + 1, 2):
                            dist = round(atoms.get_distance(p[x], p[r]), 1)
                            cross2.append(dist)
                    cross2.sort()

                    if cross1 == cross2:
                        d.append(int(j))
        tmp1 = []
        tmp2 = []
        for i in range(len(ring_tlist)):
            if i not in d:
                tmp1.append(ring_tlist[i])
                tmp2.append(ring_full[i])
        label_rings[length] = tmp1
        index_rings[length] = tmp2

    return index_rings, label_rings


def get_paths(graph: nx.Graph, index: int, ring_sizes: list) -> list:
    """Get all the paths (rings) between the index atom and its neighbor.

    Args:
        graph (nx.Graph): graph object representing zeolite framework in new larger cell
        index (int): index of the T atom to be analyzed
        ring_sizes (list): list of integers representing the sizes of rings to be analyzed

    Returns:
        paths (list): list of lists containing the paths
    """

    neighbors = list(nx.neighbors(graph, index))
    neighbor = neighbors[0]

    # next find all paths connecting index and neighbor
    paths = []
    for path in nx.all_simple_paths(graph, index, neighbor, cutoff=max(ring_sizes) - 1):
        if len(path) in ring_sizes:
            paths.append(path)  # noqa: PERF401

    return paths
