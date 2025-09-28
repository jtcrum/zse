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


def add_two_protons(
    atoms: Atoms,
    indices: int,
    oxygens: np.ndarray,
    silicons: np.ndarray,
    code: str,
    path: str | None = None,
) -> tuple[list[Atoms], list[str]]:
    """Add two protons to the structure at specified sites.

    Args:
        atoms (Atoms): The ASE Atoms object representing the structure.
        indices (int): The indices of the atoms to which the protons will be added.
        oxygens (np.ndarray): Array of indices of oxygen atoms surrounding the target atoms.
        silicons (np.ndarray): Array of indices of silicon atoms surrounding the target atoms.
        code (str): The code representing the structure type.
        path (str | None, optional): The directory path to save the modified structures.
            Defaults to None.

    Returns:
        tuple[list[Atoms], list[str]]: A list of modified Atoms objects and a list of
            location labels.
    """
    labels = site_labels(atoms, code)

    adsorbate = molecule("H")
    H_lattice = atoms + adsorbate + adsorbate

    locations = []
    traj = []
    for lidx in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[indices[0]]
        H_lattice.translate(diff)
        H_lattice.wrap()
        H_lattice.translate(-1 * diff)
        H_lattice.wrap()
        for k in range(4):
            center = H_lattice.get_center_of_mass()
            positions = atoms.get_positions()
            diff = center - positions[indices[1]]
            H_lattice.translate(diff)
            H_lattice.wrap()
            H_lattice.translate(-1 * diff)
            H_lattice.wrap()

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
