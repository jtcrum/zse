"""Tools for reading and parsing CIF files."""

import importlib.resources as pkg_resources
import os
import warnings

import numpy as np
from ase import Atoms
from ase.io import read

__all__ = ["cif_site_labels", "read_cif"]

warnings.filterwarnings("ignore")

path = ".temp_files/"
filepath = pkg_resources.files(__name__).joinpath(path)

"""
NOTE ABOUT CIF FILE FORMATS:
CIFs must include '_symmetry_Int_Taables_number' to be read by ASE.
If this is not included please edit your CIF file to include this information.
"""


def get_atom_lines(alllines: list[str]) -> tuple[int, int, list[int]]:
    """Parse the atom lines from a CIF file.

    Args:
        alllines (list[str]): List of lines from the CIF file.

    Returns:
        tuple[int, int, list[int]]: Start and end indices of atom lines and the order of
            relevant columns.
    """
    order = []
    for i, line in enumerate(alllines):
        if "_atom" in line:
            order.append(line)
            start = i + 1
    end = None
    for i, line in enumerate(alllines[start:]):
        if len(line.split()) == 0:
            end = start + i - 1
            break
    if not end:
        end = len(alllines) - 1

    new_order = []
    for i, o in enumerate(order):
        if "site_label" in o:
            new_order.append(i)
        if "site_type_symbol" in o:
            new_order.append(i)
        if "fract_x" in o:
            new_order.append(i)
        if "fract_y" in o:
            new_order.append(i)
        if "fract_z" in o:
            new_order.append(i)

    return start, end, new_order


def fix_cif(cif: str) -> tuple[Atoms, list[str]]:
    """Fix common formatting issues in a CIF file and read it with ASE.

    Args:
        cif (str): Path to the CIF file.

    Returns:
        tuple[Atoms, list[str]]: ASE Atoms object and list of lines from the
            modified CIF file.
    """
    with open(cif, "r") as f:
        alllines = f.readlines()

    for i, line in enumerate(alllines):
        if "IT_coordinate_system_code" in line:
            fields = line.split()
            alllines[i] = f"_symmetry_space_group_setting {fields[-1]} \n"

        if "_atom_site_type_symbol" in line and "_atom_site_label" in alllines[i + 1]:
            alllines[i], alllines[i + 1] = alllines[i + 1], alllines[i]

    file_name = cif.rstrip(".cif")
    temp_file = f"{filepath}/{file_name.split('/')[-1]}_temp.cif"
    with open(temp_file, "w") as f:
        f.writelines(alllines)
    atoms = read(temp_file)
    os.remove(temp_file)
    return atoms, alllines


def get_tsites(cif: str) -> tuple[list[str], list[int], list[int]]:
    """Parse the T-sites from a CIF file.

    Args:
        cif (str): Path to the CIF file.

    Returns:
        tuple[list[str], list[int], list[int]]: List of T-site labels, their multiplicities,
            and the corresponding atom indices.
    """
    tsites = []
    tpos = []
    z, alllines = fix_cif(cif)
    si = [atom.index for atom in z if atom.symbol != "O"]
    start, end, order = get_atom_lines(alllines)
    for line in alllines[start : end + 1]:
        if "Si" in line or "T" in line:
            line = line.split()
            temp_label = line[order[0]]
            if not any(str.isdigit(c) for c in temp_label):
                temp_label = line[order[1]]
            if "Si" in temp_label:
                temp_label = temp_label.replace("Si", "T")
            tsites.append(temp_label)
            pos = [float(line[order[2]]), float(line[order[3]]), float(line[order[4]])]
            tpos.append([round(num, 2) for num in pos])

    tpos = np.array(tpos)
    pos = z[si].get_scaled_positions()
    tinds = []
    for tp in tpos:
        for i, p in enumerate(pos):
            p = [round(num, 2) for num in p]
            diff = abs(tp - p)
            if sum(diff) <= 0.03:
                tinds.append(si[i])

    tmults = [tinds[i] - tinds[i - 1] for i in range(1, len(tinds))]
    tmults.append(si[-1] - tinds[-1] + 1)

    n = len(si)
    sn = sum(tmults)
    if n != sn:
        print("Something Went Wrong With T Sites")
    return tsites, tmults, tinds


def get_osites(cif: str) -> tuple[list[str], list[int], list[int]]:
    """Parse the O-sites from a CIF file.

    Args:
        cif (str): Path to the CIF file.

    Returns:
        tuple[list[str], list[int], list[int]]: List of O-site labels, their
            multiplicities, and the corresponding atom indices.
    """
    osites = []
    opos = []
    z, alllines = fix_cif(cif)
    start, end, order = get_atom_lines(alllines)
    for line in alllines[start : end + 1]:
        if "O" in line:
            line = line.split()
            temp_label = line[order[0]]
            if not any(str.isdigit(c) for c in temp_label):
                temp_label = line[order[1]]
            osites.append(temp_label)
            pos = [float(line[order[2]]), float(line[order[3]]), float(line[order[4]])]
            opos.append([round(num, 2) for num in pos])
    opos = np.array(opos)
    pos = z.get_scaled_positions()
    oinds = []

    o = [atom.index for atom in z if atom.symbol == "O"]
    o_pos = z[o].get_scaled_positions()
    for op in opos:
        for i, p in enumerate(o_pos):
            p = np.array([round(num, 2) for num in p])
            diff = abs(op - p)
            if sum(diff) <= 0.02:
                oinds.append(o[i])

    omults = [oinds[i] - oinds[i - 1] for i in range(1, len(oinds))]
    omults.append(o[-1] - oinds[-1] + 1)

    n = len(o)
    sn = sum(omults)
    if n != sn:
        print("Something Went Wrong With O Sites")
    return osites, omults, oinds


def read_cif(
    cif: str,
) -> tuple[Atoms, list[str], list[int], list[int], list[str], list[int], list[int]]:
    """Read a CIF file and extract T-site and O-site information.

    Args:
        cif (str): Path to the CIF file.

    Returns:
        tuple: Contains the ASE Atoms object, T-site labels, T-site multiplicities,
               T-site indices, O-site labels, O-site multiplicities, and O-site indices.
    """
    atoms, _all_lines = fix_cif(cif)
    ts, tm, tinds = get_tsites(cif)
    os, om, oinds = get_osites(cif)
    return atoms, ts, tm, tinds, os, om, oinds


def cif_site_labels(cif: str) -> dict[int, str]:
    """Generate a mapping of atom indices to site labels from a CIF file.

    Args:
        cif (str): Path to the CIF file.

    Returns:
        dict[int, str]: Dictionary mapping atom indices to their corresponding site labels.
    """
    _atoms, ts, tm, tinds, os, om, oinds = read_cif(cif)
    labels = {}
    for i, t in enumerate(ts):
        for j in range(tm[i]):
            labels[tinds[i] + j] = t

    for i, o in enumerate(os):
        for j in range(om[i]):
            labels[oinds[i] + j] = o

    return labels
