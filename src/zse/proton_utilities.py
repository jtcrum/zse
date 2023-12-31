from __future__ import annotations

import os
from copy import deepcopy

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import write

from zse.utilities import site_labels


def get_os_and_ts(atoms, index):
    lattice = deepcopy(atoms)
    total_oxygen = [atom.index for atom in lattice if atom.symbol == "O"]
    total_silicon = [atom.index for atom in lattice if atom.symbol == "Si"]

    oxygens = []
    for k in total_oxygen:
        distance = lattice.get_distances(index, k, mic=True)
        if distance < 2.0:
            oxygens.append(k)
    oxygens = np.array(oxygens)
    silicons = []
    for l in oxygens:
        tmp = []
        for m in total_silicon:
            distance = lattice.get_distance(l, m, mic=True)
            if distance < 2.0:
                silicons.append(m)
    silicons = np.array(silicons)

    return oxygens, silicons


def add_one_proton(atoms, index, oxygens, silicons, code, path=None):
    labels = site_labels(atoms, code)

    hydrogen = [len(atoms)]

    adsorbate = molecule("H")
    adsorbate.translate([0, 0, 0])
    H_lattice = atoms + adsorbate

    traj = []
    locations = []
    for l in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[index]
        H_lattice.translate(diff)
        H_lattice.wrap()
        H_lattice.set_distance(oxygens[l], hydrogen[0], 0.98, fix=0)
        H_lattice.set_angle(
            int(index), int(oxygens[l]), int(hydrogen[0]), 109.6, mask=None
        )
        H_lattice.set_angle(
            int(silicons[l]), int(oxygens[l]), int(hydrogen[0]), 109.6, mask=None
        )
        H_lattice.set_dihedral(
            int(index), int(oxygens[l]), int(silicons[l]), hydrogen[0], 180, mask=None
        )
        H_lattice.translate(-1 * diff)
        H_lattice.wrap()
        traj += [Atoms(H_lattice)]
        locations.append(labels[oxygens[l]])

        if path:
            os.makedirs(f"{path}/D-{labels[oxygens[l]]}", exist_ok=True)

            write(
                f"{path}/D-{labels[oxygens[l]]}/POSCAR",
                H_lattice,
                sort=False,
            )

    return traj, locations


def add_two_protons(atoms, indices, oxygens, silicons, code, path=None):
    labels = site_labels(atoms, code)

    adsorbate = molecule("H")
    H_lattice = atoms + adsorbate + adsorbate

    locations = []
    traj = []
    for l in range(4):
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
            locations.append(f"{labels[oxygens[0][l]]}-{labels[oxygens[1][k]]}")
            if path:
                os.makedirs(
                    f"{path}/D-{labels[oxygens[0][l]]}-{labels[oxygens[1][k]]}",
                    exist_ok=True,
                )

                write(
                    f"{path}/D-{labels[oxygens[0][l]]}-{labels[oxygens[1][k]]}/POSCAR",
                    H_lattice,
                    sort=False,
                )

    return traj, locations
