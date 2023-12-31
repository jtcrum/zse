from __future__ import annotations

import math
from collections import defaultdict
from copy import deepcopy

import numpy as np

from zse.collections.framework import get_framework
from zse.ring_utilities import atoms_to_graph
from zse.rings import get_unique_rings
from zse.substitute import tsub
from zse.utilities import get_tsites, site_labels


def get_pairs(code, validation=None, max_ring=12):
    z = get_framework(code)
    tinds = get_tsites(code)[2]
    c, r, _, a = get_unique_rings(z, tinds, validation=validation, max_ring=max_ring)
    lr = defaultdict(list)
    tr = defaultdict(list)

    index_list = np.arange(len(a)).reshape(-1, len(z))
    labels = site_labels(z, code)
    for p in r:
        p.insert(0, p.pop())
        tr[len(p) // 2].append(p)
        temp = []
        for x in p:
            ind = int(np.where(index_list == x)[1])
            temp.append(labels[ind])
        lr[len(p) // 2].append(temp)

    max_ring *= 2
    traj = []
    alltlist = []
    allpairlist = []
    ring_list = []
    for r in sorted(tr):
        for nn in range(2, math.floor(r / 2) + 1):
            for q, trs in enumerate(tr[r]):
                pair_list = []
                tp = trs + trs
                lp = lr[r][q] + lr[r][q]
                z = get_framework(code)
                repeat = atoms_to_graph(z, tp[0], max_ring)[2]
                z2 = z.repeat(repeat)
                for i in range(1, len(trs) - nn, 2):
                    pair_inner = lp[i - 1 : i + 2 * nn + 2]
                    pair_id = "_".join(pair_inner)
                    pair_inner.reverse()
                    r_pair_id = "_".join(pair_inner)

                    tinds = [tp[i], tp[i + 2 * nn]]

                    flag = False
                    if pair_id in pair_list or r_pair_id in pair_list:
                        flag = True

                    if not flag:
                        zl = len(z2)
                        ozl = len(z)
                        indices = np.arange(zl)
                        indices = indices.reshape(np.prod(repeat), ozl)
                        newinds = [int(np.where(indices == ti)[1]) for ti in tinds]
                        z3 = deepcopy(z)
                        z3 = tsub(z3, newinds, "Al")
                        traj += [z3]
                        pair_list.append(pair_id)
                        alltlist.append(tinds)
                        t1label = [lp[i - 1], lp[i], lp[i + 1]]
                        t2label = [
                            lp[i + 2 * nn - 1],
                            lp[i + 2 * nn],
                            lp[i + 2 * nn + 1],
                        ]

                        allpairlist.append([t1label, t2label])
                        ring_list.append(f"{r}-MR {nn}NN")
                        if nn == r / 2:
                            pair_inner = lp[i + 2 * nn - 1 : 2 * (i + 2 * nn) + 1]
                            pair_id = "_".join(pair_inner)
                            pair_inner.reverse()
                            r_pair_id = "_".join(pair_inner)
                            pair_list.append(pair_id)
    pairs = [f"{r} | {c}" for r, c in zip(ring_list, allpairlist)]
    return pairs, traj
