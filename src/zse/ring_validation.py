"""
This module contains all the various ring validation techniques implemented by
ZSE. This is a work in progress, and more methods will be added.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy

import numpy as np

from zse.ring_utilities import is_valid, shortest_valid_path
from zse.utilities import scale_cell


def sp(G, paths):
    """
    Custum method for determing valid rings by ensuring for paths of 8 T-sites
    or larger, the current path is the shortest possible way to connect each
    node long the path. Alternate paths must be shorter than the current path
    to be considered invalid.

    7 T-site paths and smaller are handled slightly
    different where if the alrternate path is equal in length to the current
    path it is considered non valid.
    """

    import networkx as nx

    valid_paths = []
    for p in paths:
        l = len(p)
        flag = True
        for j in range(1, l - 3, 2):
            if l < 16:
                for k in range(j + 4, l, 2):
                    shortest_path = nx.shortest_path(G, p[j], p[k])
                    for r in shortest_path:
                        if r not in p:
                            lengths = [k - j + 1, l - k + j + 1]
                            if len(shortest_path) < min(lengths):
                                flag = False
                                break
            else:
                for k in range(j + 2, l, 2):
                    shortest_path = nx.shortest_path(G, p[j], p[k])
                    for r in shortest_path:
                        if r not in p:
                            lengths = [k - j + 1, l - k + j + 1]
                            if len(shortest_path) < max(lengths):
                                flag = False
                                break
        if flag:
            valid_paths.append(p)

    return valid_paths


def sphere(atoms, paths, cutoff):
    """
    Method for determining valid rings by ensuring that non ring atoms are not
    within a cutoff distance of the center of mass of the path.
    """
    from ase.geometry import get_distances

    valid_paths = []
    for p in paths:
        flag = True
        if len(p) > 10:
            ring_atoms = atoms[p]
            com = ring_atoms.get_center_of_mass()
            positions = atoms.get_positions()
            distances = get_distances(com, positions)[1][0]

            for i, d in enumerate(distances):
                if d < cutoff and i not in p:
                    flag = False
                    break
        if flag:
            valid_paths.append(p)

    return valid_paths


def cross_distance(atoms, paths):
    """
    This validation method uses cross Si-Si distances of the path to determine
    if the path is a valid ring or not. The cutoffs for the cross Si-Si
    distances were determined by analyzing many valid and invalid rings.
    This method is biased by human interpretation of what constitutes a ring.
    """
    import math

    atoms = scale_cell(atoms)

    delete = []
    for j, r in enumerate(paths):
        n = len(r) // 2
        cutoff = n / 4
        cutoff = max(cutoff, 2)
        if n % 2 == 0 and n > 5:
            distances = []
            inner_flag = False
            for x in range(1, n, 2):
                dist = atoms.get_distance(r[x], r[x + n], mic=True)
                distances.append(dist)
                if dist < n - cutoff:
                    delete.append(j)
                    inner_flag = True
                    break
            if inner_flag is False:
                outer_flag = any(d > n - math.floor(n / 6) for d in distances)
                if not outer_flag:
                    delete.append(j)
        if n % 2 != 0 and n > 5:
            r2 = deepcopy(r)
            r2.append(r[:2])
            for x in range(1, n, 2):
                dist1 = atoms.get_distance(r2[x], r2[x + n - 1], mic=True)
                dist2 = atoms.get_distance(r2[x], r2[x + n + 1], mic=True)
                if dist1 < n - cutoff or dist2 < n - cutoff:
                    delete.append(j)
                    break

    tmp_paths = deepcopy(paths)
    paths = [r for j, r in enumerate(tmp_paths) if j not in delete]
    return paths


def goetzke(G, index, cutoff):
    """
    Method to find all cycles that cannot be decomposed into smaller cycles
    via shortcuts. This rule and algorithm was presented by:
    K. Goetzke and H.-J. Klein (https://doi.org/10.1016/0022-3093(91)90145-V)
    Slight modifications to the algorithm have been implemented to take
    advantage of symmetry.
    """
    import networkx as nx

    sp = dict(nx.all_pairs_shortest_path_length(G, cutoff / 2))
    paths = []
    for size in np.arange(6, cutoff + 1, 2):
        for i in sp[index]:
            x = sp[index][i]
            if x == (size / 2):
                if path := make_path(sp, index, i, size):
                    paths += path
    return paths


def get_left(cl, s1, count, sp):
    left_options = []
    for j1 in sp[cl]:
        q = sp[cl][j1]
        try:
            if q == 1 and sp[j1][s1] == count:
                left_options.append(j1)
        except Exception:
            pass
    return left_options


def get_right(cr, s2, lo, count, size, sp):
    right_options = []
    for l in lo:
        temp = []
        flag = False
        for j2 in sp[cr]:
            try:
                q = sp[cr][j2]
                if q == 1 and sp[j2][s2] == count and sp[l][j2] == size / 2:
                    temp.append(j2)
                    flag = True
            except Exception:
                pass
        if flag:
            right_options.append(temp)
        else:
            right_options.append("x")
    return right_options


def make_path(sp, s1, s2, size):
    count = size / 2 - 1
    left_overall = [[s2]]
    right_overall = [[s1]]
    while count > 0:
        left_temp = []
        right_temp = []
        for left, right in zip(left_overall, right_overall):
            cl = left[-1]
            cr = right[-1]
            left_options = get_left(cl, s1, count, sp)
            right_options = get_right(cr, s2, left_options, count, size, sp)

            for l, ro in zip(left_options, right_options):
                for r in ro:
                    if r != "x":
                        left_temp.append(left + [l])
                        right_temp.append(right + [r])
        left_overall = []
        right_overall = []
        for l, r in zip(left_temp, right_temp):
            l1 = l[-1]
            r1 = r[-1]
            flag = True

            for l2, r2 in zip(left_overall, right_overall):
                if size / 2 % 2 == 0:
                    if l1 == r2[-1] and r1 == l2[-1]:
                        flag = False
                elif l1 == r2[-2] and r1 == l2[-2]:
                    flag = False
            if flag:
                left_overall.append(l)
                right_overall.append(r)
        count -= 1
    paths = []
    for l, r in zip(left_overall, right_overall):
        path = r + l
        if len(path) == size:
            paths.append(path)
    return paths


def sastre(G, paths, index_symbol):
    """
    This method returns only the shortest paths rings based on the rule
    presented by Sastre and Corma.
    Sastre, G; Corma, A. (DOI: 10.1021/jp8100128)
    A valid ring is a cycle that cannot be decomposed into smaller rings, and
    is considered a vertex symbol ring for at least one pair of nearest
    neighbor oxygens in the cycle.
    Results found with this method match results from the Sastre & Corma paper.

    """

    start = 0 if index_symbol == "O" else 1
    valid_paths = []
    for p in paths:
        p2 = p + p[:4]
        l = len(p)
        for j in range(start, len(p), 2):
            flag = False
            path, length = shortest_valid_path(G, p2[j], p2[j + 2], p2[j + 1], l)
            if length == l:
                flag = True
                break
        if flag:
            valid_paths.append(p)
    return valid_paths


def crum(G, paths, index_symbol):
    """
    Method to remove composite stacked rings (i.e. 8-MRs in the CHA D6R) or
    14-MRs in the AFI framework.
    """
    import networkx as nx

    start = 1 if index_symbol == "O" else 0
    valid_paths = []
    for path in paths:
        FLAG = False
        path2 = path + path
        l = len(path)
        if l > 8 and (l / 2) % 2 == 0:
            for j in range(start, l // 2 - 1, 2):
                for k in [int(j + l / 2 - 2)]:
                    p1 = path2[j : k + 1]
                    p2 = path2[k : l + j + 1]
                    l1 = len(p1)
                    l2 = len(p2)
                    if l1 < l2:
                        G2 = deepcopy(G)
                        for x in p1[1:-1]:
                            G2.remove_node(x)
                        sp = nx.shortest_path(G2, p1[0], p1[-1])
                        con1 = p2 + sp[1:-1]
                        if len(con1) < l:
                            FLAG = True
                        if len(con1) == l:
                            FLAG, _ = is_valid(G, con1)
                    if l2 < l1:
                        G2 = deepcopy(G)
                        for x in p2[1:-1]:
                            G2.remove_node(x)
                        sp = nx.shortest_path(G2, p2[0], p2[-1])
                        con1 = p1 + sp[1:-1]
                        if len(con1) < l:
                            FLAG = True
                        if len(con1) == l:
                            FLAG, _ = is_valid(G, con1)
                    if l1 == l2:
                        G2 = deepcopy(G)
                        sp = nx.shortest_path(G2, p1[0], p1[-1])
                        if len(sp) < l1:
                            FLAG = True
                    if FLAG:
                        break
                if FLAG:
                    break
        if not FLAG:
            valid_paths.append(path)
    return valid_paths


def vertex(paths):
    oxygens = []
    v_paths = defaultdict(list)
    for p in paths:
        oxygens.extend((p[1], p[-1]))
        sites = [p[1], p[-1]]
        sites.sort()
        v_paths[f"{sites[0]}-{sites[1]}"].append(p)
    oxygens = np.unique(oxygens)
    valid_paths = []
    for v in sorted(v_paths):
        paths = v_paths[v]
        l = len(paths[0])
        valid_paths.extend(p for p in paths if len(p) == l)
    return valid_paths
