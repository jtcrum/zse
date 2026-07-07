"""Utilities for adding protons to structures."""

import os

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import write

from zse.utilities import site_labels

__all__ = ["add_one_proton", "add_two_protons", "get_os_and_ts"]


def get_os_and_ts(atoms: Atoms, index: int) -> tuple[np.ndarray, np.ndarray]:
    """Get the oxygen and silicon atoms surrounding a given index.

    Args:
        atoms (Atoms): The ASE Atoms object representing the structure.
        index (int): The index of the atom for which to find surrounding oxygens and silicons.

    Returns:
        tuple[np.ndarray, np.ndarray]: Two arrays containing the indices of the surrounding
            oxygen and silicon atoms, respectively.
    """
    lattice = atoms.copy()
    total_oxygen = [atom.index for atom in lattice if atom.symbol == "O"]
    total_silicon = [atom.index for atom in lattice if atom.symbol == "Si"]

    oxygens = []
    for k in total_oxygen:
        distance = lattice.get_distances(index, k, mic=True)
        if distance < 2.0:
            oxygens.append(k)
    oxygens = np.array(oxygens)
    silicons = []
    for lidx in oxygens:
        for midx in total_silicon:
            distance = lattice.get_distance(lidx, midx, mic=True)
            if distance < 2.0:
                silicons.append(midx)
    silicons = np.array(silicons)

    return oxygens, silicons


def add_one_proton(
    atoms: Atoms,
    index: int,
    oxygens: np.ndarray,
    silicons: np.ndarray,
    code: str,
    path: str | None = None,
) -> tuple[list[Atoms], list[str]]:
    """Add a single proton to the structure at specified sites.

    Args:
        atoms (Atoms): The ASE Atoms object representing the structure.
        index (int): The index of the atom to which the proton will be added.
        oxygens (np.ndarray): Array of indices of oxygen atoms surrounding the target atom.
        silicons (np.ndarray): Array of indices of silicon atoms surrounding the target atom.
        code (str): The code representing the structure type.
        path (str | None, optional): The directory path to save the modified structures.
            Defaults to None.

    Returns:
        tuple[list[Atoms], list[str]]: A list of modified Atoms objects and a list of
            location labels.
    """
    labels = site_labels(atoms, code)

    hydrogen = [len(atoms)]

    adsorbate = molecule("H")
    adsorbate.translate([0, 0, 0])
    H_lattice = atoms + adsorbate

    traj = []
    locations = []
    for lidx in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[index]
        H_lattice.translate(diff)
        H_lattice.wrap()
        H_lattice.set_distance(oxygens[lidx], hydrogen[0], 0.98, fix=0)
        H_lattice.set_angle(int(index), int(oxygens[lidx]), int(hydrogen[0]), 109.6, mask=None)
        H_lattice.set_angle(
            int(silicons[lidx]), int(oxygens[lidx]), int(hydrogen[0]), 109.6, mask=None
        )
        H_lattice.set_dihedral(
            int(index), int(oxygens[lidx]), int(silicons[lidx]), hydrogen[0], 180, mask=None
        )
        H_lattice.translate(-1 * diff)
        H_lattice.wrap()
        traj += [Atoms(H_lattice)]
        locations.append(labels[oxygens[lidx]])

        if path:
            os.makedirs(f"{path}/D-{labels[oxygens[lidx]]}", exist_ok=True)

            write(f"{path}/D-{labels[oxygens[lidx]]}/POSCAR", H_lattice, sort=False)

    return traj, locations


def _place_proton(
    H_lattice: Atoms, t_index: int, o_index: int, si_index: int, h_index: int, dihedral: float = 180
) -> None:
    """Bond a single proton (already present in H_lattice at h_index) to o_index,
    in a tetrahedral arrangement about the T-O-Si linkage. Modifies H_lattice in place.
    """
    center = H_lattice.get_center_of_mass()
    positions = H_lattice.get_positions()
    diff = center - positions[t_index]
    H_lattice.translate(diff)
    H_lattice.wrap()
    H_lattice.set_distance(int(o_index), h_index, 0.98, fix=0)
    H_lattice.set_angle(int(t_index), int(o_index), h_index, 109.6, mask=None)
    H_lattice.set_angle(int(si_index), int(o_index), h_index, 109.6, mask=None)
    H_lattice.set_dihedral(int(t_index), int(o_index), int(si_index), h_index, dihedral, mask=None)
    H_lattice.translate(-1 * diff)
    H_lattice.wrap()


def _clashes(
    H_lattice: Atoms, h_index: int, bonded_o_index: int, min_distance: float = 1.0
) -> bool:
    """Check whether the proton at h_index sits too close to any other atom
    (other than the oxygen it's bonded to, which is intentionally 0.98 A away).
    """
    others = [i for i in range(len(H_lattice)) if i not in (h_index, int(bonded_o_index))]
    distances = H_lattice.get_distances(h_index, others, mic=True)
    return bool(np.any(distances < min_distance))


def add_two_protons(
    atoms: Atoms,
    indices: list[int],
    oxygens: np.ndarray,
    silicons: np.ndarray,
    code: str,
    path: str | None = None,
) -> tuple[list[Atoms], list[str]]:
    """Add two protons to the structure, enumerating every combination of the
    oxygens surrounding each of the two T-sites (e.g. 4x4=16 structures for
    two 4-coordinate T-sites). Each proton is placed with the same tetrahedral
    T-O-H/Si-O-H geometry used by 'add_one_proton'; combinations where a
    proton would clash with another atom (including the other new proton) are
    resolved by flipping that proton's dihedral to the other staggered position.

    Args:
        atoms (Atoms): The ASE Atoms object representing the structure.
        indices (list[int]): The indices of the two T-sites to protonate.
        oxygens (np.ndarray): The two arrays of oxygen indices surrounding each T-site.
        silicons (np.ndarray): The two arrays of Si indices bonded to those oxygens.
        code (str): The code representing the structure type.
        path (str | None, optional): The directory path to save the modified structures.
            Defaults to None.

    Returns:
        tuple[list[Atoms], list[str]]: A list of modified Atoms objects and a list of
            location labels.
    """
    labels = site_labels(atoms, code)

    adsorbate = molecule("H")
    base_lattice = atoms + adsorbate + adsorbate
    h1, h2 = len(atoms), len(atoms) + 1

    traj = []
    locations = []
    for lidx in range(4):
        for k in range(4):
            H_lattice = base_lattice.copy()
            _place_proton(H_lattice, indices[0], oxygens[0][lidx], silicons[0][lidx], h1)
            _place_proton(H_lattice, indices[1], oxygens[1][k], silicons[1][k], h2)

            if _clashes(H_lattice, h1, oxygens[0][lidx]) or _clashes(H_lattice, h2, oxygens[1][k]):
                # retry the second proton at the other staggered position
                _place_proton(
                    H_lattice, indices[1], oxygens[1][k], silicons[1][k], h2, dihedral=0
                )

            traj += [Atoms(H_lattice)]
            locations.append(f"{labels[oxygens[0][lidx]]}-{labels[oxygens[1][k]]}")
            if path:
                os.makedirs(
                    f"{path}/D-{labels[oxygens[0][lidx]]}-{labels[oxygens[1][k]]}",
                    exist_ok=True,
                )

                write(
                    f"{path}/D-{labels[oxygens[0][lidx]]}-{labels[oxygens[1][k]]}/POSCAR",
                    H_lattice,
                    sort=False,
                )

    return traj, locations
