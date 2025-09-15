__all__ = ["add_cation", "count_rings"]

import os

import numpy as np
from ase import Atoms
from ase.io import write


def count_rings(paths):
    _class = [int(len(p) / 2) for p in paths]
    paths = [x for _, x in sorted(zip(_class, paths, strict=False), reverse=True)]
    _class.sort(reverse=True)

    class_count = []
    for i in range(len(_class)):
        if i == 0:
            class_count.append(1)
        else:
            counter = 1
            for j in range(i):
                if _class[i] == _class[j]:
                    counter += 1
            class_count.append(counter)

    return _class, class_count, paths


def add_cation(
    atoms,
    large_atoms,
    radii,
    index,
    symbol,
    paths,
    included_rings,
    class_count,
    path=None,
    bvect=None,
):
    traj = []
    locations = []
    for i in range(len(paths)):
        p = paths[i]
        if len(p) in included_rings:
            trans = atoms[index].position - large_atoms[index].position
            positions = large_atoms[p].positions
            co = sum(positions) / len(positions)
            co += trans
            vector = co - atoms[index].position
            vector = np.array(vector)
            vhat = vector / np.linalg.norm(vector)

            bond_length = bvect or radii[symbol] + 2 / radii[symbol]

            new_p = [atoms[index].position + bond_length * vhat]

            adsorbate = Atoms(symbol)
            adsorbate.set_positions(new_p)
            c_atoms = atoms + adsorbate
            c_atoms.wrap()
            traj += [c_atoms]
            locations.append(f"{int(len(p) / 2)!s}MR")
            # write POSCAR for each structure

            if path:
                os.makedirs(
                    f"{path}/D-{int(len(p) / 2)!s}MR-{class_count[i]!s}",
                    exist_ok=True,
                )

                write(
                    f"{path}/D-{int(len(p) / 2)!s}MR-{class_count[i]!s}/POSCAR",
                    c_atoms,
                    sort=True,
                )
    return traj, locations
