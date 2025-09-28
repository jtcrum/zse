"""Utilities for adding cations to zeolite frameworks and counting rings."""

import os
from typing import Iterable

import numpy as np
from ase import Atoms
from ase.io import write

__all__ = ["add_cation", "count_rings"]


def count_rings(paths: Iterable[str]) -> tuple[list[int], list[int], list[str]]:
    """Count the number of rings of each size in a list of paths.

    Args:
        paths (Iterable[str]): Iterable of paths representing rings.

    Returns:
        list[int]: A list of ring sizes
        list[int]: A list of counts for each ring size
        list[str]: The original list of paths.
    """
    _class = [int(len(p) / 2) for p in paths]
    paths = [x for _, x in sorted(zip(_class, paths, strict=True), reverse=True)]
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
    atoms: Atoms,
    large_atoms: Atoms,
    radii: dict[str, float],
    index: int,
    symbol: str,
    paths: Iterable[Iterable[int]],
    included_rings: set[int],
    class_count: list[int],
    path: str | None = None,
    bvect: np.ndarray | None = None,
) -> tuple[list[Atoms], list[str]]:
    """Add a cation to an atoms object at specified ring locations.

    Args:
        atoms (Atoms): The ASE Atoms object of the zeolite framework.
        large_atoms (Atoms): The ASE Atoms object of the larger repeated framework.
        radii (dict[str, float]): A dictionary mapping element symbols to their atomic radii.
        index (int): The index of the atom in the Atoms object where the cation will be added.
        symbol (str): The symbol of the cation to be added (e.g., 'Na', 'K').
        paths (Iterable[Iterable[int]]): A list of paths representing rings in the framework.
        included_rings (set[int]): A set of ring sizes to include for cation placement.
        class_count (list[int]): A list of counts for each ring size class.
        path (str | None, optional): The directory path to save POSCAR files. Defaults to None.
        bvect (np.ndarray | None, optional): A bond vector for cation placement. Defaults to None.

    Returns:
        list[Atoms]: a list of Atoms objects with the cation added
        list[str]: a list of location labels for each cation placement.
    """
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
