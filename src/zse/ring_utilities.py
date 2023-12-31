"""
This module contains utilities to be used by the rings.py module.
"""
from __future__ import annotations

from collections import defaultdict as dd
from copy import deepcopy
from itertools import permutations as perm

import networkx as nx
import numpy as np

from zse.utilities import center


def atoms_to_graph(atoms, index, max_ring):
    """
    Helper function to repeat a unit cell enough times to capture the largest
    possible ring, and turn the new larger cell into a graph object.

    RETURNS:
    G = graph object representing zeolite framework in new larger cell
    large_atoms = ASE atoms object of the new larger cell framework
    repeat = array showing the number of times the cell was repeated: [x,y,z]
    """

    # repeat cell, center the cell, and wrap the atoms back into the cell
    cell = atoms.cell.cellpar()[:3]
    repeat = []
    for c in cell:
        if c / 2 < max_ring / 2 + 5:
            l = c
            re = 1
            while l / 2 < max_ring / 2 + 5:
                re += 1
                l = c * re

            repeat.append(re)
        else:
            repeat.append(1)
    large_atoms = deepcopy(atoms)
    large_atoms = large_atoms.repeat(repeat)
    center = large_atoms.get_center_of_mass()
    trans = center - large_atoms.positions[index]
    large_atoms.translate(trans)
    large_atoms.wrap()

    # remove atoms that won't contribute to wrings
    from ase.geometry import get_distances

    cell = large_atoms.get_cell()
    p1 = large_atoms[index].position
    positions = large_atoms.get_positions()
    distances = get_distances(p1, positions)[1][0]

    delete = [i for i, l in enumerate(distances) if l > max_ring / 2 + 5]
    inds = [atom.index for atom in large_atoms]
    large_atoms.set_tags(inds)
    atoms = deepcopy(large_atoms)
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
    # now we make the graph
    import networkx as nx

    G = nx.from_numpy_array(matrix)
    # G.remove_nodes_from(delete)
    return G, large_atoms, repeat


def remove_dups(paths):
    """
    This is a helper function for get_orings and get_trings.
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


def paths_to_atoms(atoms, paths):
    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms if atom.index not in keepers]
    tmp_atoms = deepcopy(atoms)
    del tmp_atoms[d]

    return tmp_atoms


def remove_geometric_dups(atoms, paths):
    unique_paths = []
    unique_distances = []
    for p in paths:
        atoms, _ = center(atoms, p[0])
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


def get_vertices(G, index):
    vertices = []
    neighbors = list(nx.neighbors(G, index))
    l = len(neighbors)
    for j in range(l - 1):
        vertices.extend([neighbors[j], neighbors[k]] for k in range(j + 1, l))
    return vertices


def shortest_valid_path(G, o1, o2, index, l):
    import networkx as nx

    G2 = deepcopy(G)
    G2.remove_node(index)
    flag = True
    l2 = 6
    p = []
    while flag:
        paths = nx.all_simple_paths(G2, o1, o2, l2 - 1)
        for p in paths:
            p.append(index)
            if len(p) == l2:
                flag, j = is_valid(G, p)
            if not flag:
                break
        if flag:
            l2 += 2
        if l2 > l:
            p = [1]
            break
    return p, len(p)


def is_valid(G, path):
    import networkx as nx

    l = len(path)
    flag = False
    for j in range(1, l - 1, 2):
        for k in range(j + 2, l, 2):
            node = path[j]
            sp = nx.shortest_path(G, node, path[k])

            if len(sp) < k - j + 1 and len(sp) < l - (k - j) + 1:
                flag = True
                break
        if flag:
            break
    return flag, j


def all_paths(G, o1, o2, index, l):
    all_paths = []
    G2 = deepcopy(G)
    G2.remove_node(index)
    paths = nx.all_simple_paths(G2, o1, o2, l)
    for path in paths:
        path.append(index)
        if len(path) == l:
            flag, j = is_valid(G, path)
            if not flag and path not in all_paths:
                all_paths.append(path)
    return all_paths


def vertex_order(r):
    oxygens = []
    o_pair_paths = dd(list)
    o_pair_sizes = dd(int)
    o_pair_counts = dd(lambda: 0)
    o_pair_weights = dd(lambda: 0)

    for path in r:
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
            except Exception:
                pass
        weights.append(w)

    zipped_lists = zip(weights, perms)
    sp = sorted(zipped_lists, reverse=True)
    tuples = zip(*sp)
    weights, perms = (list(tuple) for tuple in tuples)

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
            new_r.extend(iter(o_pair_paths[k]))
        except Exception:
            pass

    ordered_v = []
    for c, s in zip(counts, sizes):
        if c == 0:
            ordered_v.append("*")
        elif c == 1:
            ordered_v.append(f"{int(s / 2)}")
        else:
            ordered_v.append(f"{int(s / 2)}_{c}")
    ordered_v = "•".join(ordered_v)
    c = [len(x) // 2 for x in new_r]

    return ordered_v, new_r


def remove_labeled_dups(index_paths, label_paths, ring_sizes, atoms):
    # first make dictionaries
    label_rings = {}
    index_rings = {}
    for i, r in enumerate(label_paths):
        length = len(r) // 2
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
                if st2 in f"{st1} {st1}" or st2_2 in f"{st1} {st1}":
                    p = ring_full[i]
                    atoms, trans = center(atoms, p[0])
                    cross1 = []
                    for x in range(1, len(p) + 1, 2):
                        for r in range(x + 2, len(p) + 1, 2):
                            dist = round(atoms.get_distance(p[x], p[r]), 1)
                            cross1.append(dist)
                    cross1.sort()

                    p = ring_full[j]
                    atoms, trans = center(atoms, p[0])
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


def get_paths(G, index, ring_sizes):
    """
    Get all the paths (rings) between the index atom and its neighbor
    """

    # first find the neighbor
    import networkx as nx

    neighbors = list(nx.neighbors(G, index))
    neighbor = neighbors[0]

    return [
        path
        for path in nx.all_simple_paths(G, index, neighbor, cutoff=max(ring_sizes) - 1)
        if len(path) in ring_sizes
    ]


def dict_to_atoms(index_paths, atoms):
    trajectories = {}
    com = atoms.get_center_of_mass()
    for length in index_paths.keys():
        paths = index_paths[length]
        tmp_traj = []
        for p in paths:
            tmp_atoms = paths_to_atoms(atoms, [p])
            trans = com - tmp_atoms[0].position
            tmp_atoms.translate(trans)
            tmp_atoms.wrap()
            tmp_traj += [tmp_atoms]
        trajectories[length] = tmp_traj

    return trajectories
