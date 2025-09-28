"""This module contains general utilities
to be used by other functions in ZSE.
"""

import numpy as np
from ase import Atoms
from ase.geometry import get_distances

from zse.collections.framework import framework_db, get_osites_db, get_tsites_db

__all__ = ["center", "get_osites", "get_tsites", "scale_cell", "site_labels"]

SILICON_DISTANCE = 3.1


def center(atoms: Atoms, index: int) -> tuple[Atoms, np.ndarray]:
    """Center an atoms object around an index.

    Args:
        atoms (Atoms): Atoms object to be centered.
        index (int): Index of the atom to center around.

    Returns:
        tuple[Atoms, np.ndarray]: Centered atoms object and translation vector.
    """
    com = atoms.get_center_of_mass()
    trans = com - atoms[index].position
    atoms.translate(trans)
    atoms.wrap()
    return atoms, trans


def get_osites(code: str) -> tuple[list[str], list[int], list[int]]:
    """Get the oxygen sites for a given zeolite framework code.

    Args:
        code (str): Zeolite framework code.

    Returns:
        tuple[list[str], list[int], list[int]]: A tuple containing a list of oxygen site labels,
        a list of their multiplicities, and a list of the first oxygen atom indices for each site.
    """
    z = framework_db(code)
    osites, omult = get_osites_db(code)
    oinds = [atom.index for atom in z if atom.symbol == "O"]
    index = 0
    first_os = []
    for m in omult:
        first_os.append(oinds[index])
        index += m
    return osites, omult, first_os


def get_tsites(code: str) -> tuple[list[str], list[int], list[int]]:
    """Get the T sites for a given zeolite framework code.

    Args:
        code (str): Zeolite framework code.

    Returns:
        list: T-site labels
        list: T-site multiplicities
        list: the first T atom indices for each site
    """
    z = framework_db(code)
    tsites, tmult = get_tsites_db(code)
    tinds = [atom.index for atom in z if atom.symbol != "O"]
    index = 0
    first_ts = []
    for m in tmult:
        first_ts.append(tinds[index])
        index += m
    return tsites, tmult, first_ts


def label_osites(atoms: Atoms, code: str) -> dict[int, str]:
    """Label the oxygen sites in an atoms object based on a zeolite framework code.

    Args:
        atoms (Atoms): Atoms object containing the zeolite framework.
        code (str): Zeolite framework code.

    Returns:
        dict[int, str]: A dictionary mapping atom indices to their corresponding oxygen site labels.
    """
    z = framework_db(code)
    osites, omult, _first = get_osites(code)

    zcell = z.cell.cellpar()[:3]
    acell = atoms.cell.cellpar()[:3]
    repeat = []
    for zc, ac in zip(zcell, acell, strict=True):
        repeat.append(int(round(ac / zc)))  # noqa: RUF046

    z = z.repeat(repeat)
    oinds = [atom.index for atom in z if atom.symbol == "O"]

    rp = np.prod(repeat)
    dict_ = {}
    j = 0
    for _i in range(rp):
        for s, t in enumerate(osites):
            for _q in range(omult[s]):
                dict_[oinds[j]] = t
                j += 1

    return dict_


def label_tsites(atoms: Atoms, code: str) -> dict[int, str]:
    """Label the T sites in an atoms object based on a zeolite framework code.

    Args:
        atoms (Atoms): Atoms object containing the zeolite framework.
        code (str): Zeolite framework code.

    Returns:
        dict[int, str]: A dictionary mapping atom indices to their corresponding T site labels.
    """
    z = framework_db(code)
    tsites, tmult, _first = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol != "O"]

    zcell = z.cell.cellpar()[:3]
    acell = atoms.cell.cellpar()[:3]
    repeat = []
    for zc, ac in zip(zcell, acell, strict=True):
        repeat.append(int(round(ac / zc)))  # noqa: RUF046
    z = z.repeat(repeat)
    tinds = [atom.index for atom in z if atom.symbol != "O"]

    rp = np.prod(repeat)
    dict_ = {}
    j = 0
    for _ in range(rp):
        for s, t in enumerate(tsites):
            for _ in range(tmult[s]):
                dict_[tinds[j]] = t
                j += 1

    return dict_


def scale_cell(atoms: Atoms) -> Atoms:
    """Scale the cell of an atoms object to have a Si-Si distance of 3.1 Angstroms.

    Args:
        atoms (Atoms): Atoms object containing the zeolite framework.

    Returns:
        Atoms: Scaled atoms object.
    """
    diff = 1
    mult = 1
    si = [atom.index for atom in atoms if atom.symbol == "Si"]
    zsi = atoms[si]
    while diff > 0.01:
        cell = atoms.cell.cellpar()
        for i in range(3):
            cell[i] = cell[i] * mult
        zsi.set_cell(cell, scale_atoms=True)
        ncell = zsi.get_cell()
        positions = zsi.get_positions()
        distances = get_distances(positions, cell=ncell, pbc=[1, 1, 1])[1]
        temp = []
        for line in distances:
            masked_line = np.ma.masked_equal(line, 0.0, copy=False)
            temp.append(masked_line.min())
        silm = np.average(temp)
        diff = abs(silm - SILICON_DISTANCE)
        mult = SILICON_DISTANCE / silm
    atoms.set_cell(cell, scale_atoms=True)

    return atoms


def site_labels(atoms: Atoms, code: str) -> dict[int, str]:
    """Get the atom site labels (as defined by the IZA) for
    your atoms object. Be sure to remove any adsorbates from your zeolite
    framework before using this function or else it won't work. This function
    will work with T sites that have been exchanged for Al (or any atom).

    Args:
        atoms (Atoms): Atoms object containing a zeolite you want labels for.
        code (str): The zeolite framework code of your atoms object (i.e. 'CHA').

    Returns:
        dict[int, str]: A dictionary mapping atom indices to their corresponding site labels.
    """

    tdict = label_tsites(atoms, code)
    odict = label_osites(atoms, code)
    all_labels = {**tdict, **odict}

    z = framework_db(code)
    zcell = z.cell.cellpar()[:3]
    acell = atoms.cell.cellpar()[:3]
    repeat = []
    for zc, ac in zip(zcell, acell, strict=True):
        repeat.append(int(round(ac / zc)))  # noqa: RUF046
    z = z.repeat(repeat)

    zo_inds = [atom.index for atom in z if atom.symbol == "O"]
    zt_inds = [atom.index for atom in z if atom.symbol != "O"]

    # z.set_tags(z_inds)
    poszo = z[zo_inds].get_scaled_positions()
    poszt = z[zt_inds].get_scaled_positions()

    scaledp = atoms.get_scaled_positions()

    dict_ = {}
    for a in atoms:
        pa = scaledp[a.index]
        sym = a.symbol
        if sym == "O":
            diffp = poszo - pa
            mags = [np.linalg.norm(d) for d in diffp]
            ind = mags.index(min(mags))
            ind = zo_inds[ind]
            label = all_labels[ind]

        if sym != "O":
            diffp = poszt - pa
            mags = [np.linalg.norm(d) for d in diffp]
            ind = mags.index(min(mags))
            ind = zt_inds[ind]
            label = all_labels[ind]
        dict_[a.index] = label

    return dict_
