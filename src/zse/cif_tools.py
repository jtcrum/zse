from __future__ import annotations

import os

import numpy as np
import pkg_resources
from ase.io import read

path = ".temp_files/"
filepath = pkg_resources.resource_filename(__name__, path)

"""
NOTE ABOUT CIF FILE FORMATS:
CIFs must include '_symmetry_Int_Taables_number' to be read by ASE.
If this is not included please edit your CIF file to include this information.
"""


def get_atom_lines(alllines):
    order = []
    for i, line in enumerate(alllines):
        if "_atom" in line:
            order.append(line)
            start = i + 1
    end = (
        next(
            (
                start + i - 1
                for i, line in enumerate(alllines[start:])
                if len(line.split()) == 0
            ),
            None,
        )
        or len(alllines) - 1
    )

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


def fix_cif(cif):
    with open(cif) as f:
        alllines = f.readlines()
    for i, line in enumerate(alllines):
        if "IT_coordinate_system_code" in line:
            fields = line.split()
            alllines[i] = f"_symmetry_space_group_setting {fields[-1]} \n"

        if "_atom_site_type_symbol" in line and "_atom_site_label" in alllines[i + 1]:
            alllines[i], alllines[i + 1] = alllines[i + 1], alllines[i]

    file_name = cif.rstrip(".cif")
    temp_file = "{0}/{1}_temp.cif".format(filepath, file_name.split("/")[-1])
    with open(temp_file, "w") as f:
        f.writelines(alllines)
    atoms = read(temp_file)
    os.remove(temp_file)
    return atoms, alllines


def get_tsites(cif):
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
    t_class = []
    for tp in tpos:
        for i, p in enumerate(pos):
            p = [round(num, 2) for num in p]
            diff = abs(tp - p)
            if sum(diff) <= 0.03:
                tinds.append(si[i])

    tmults = [tinds[i] - tinds[i - 1] for i in range(1, len(tsites))]
    tmults.append(si[-1] - tinds[-1] + 1)

    n = len(si)
    sn = sum(tmults)
    if n != sn:
        raise ValueError("Something Went Wrong With T Sites")
    return tsites, tmults, tinds


def get_osites(cif):
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

    omults = [oinds[i] - oinds[i - 1] for i in range(1, len(osites))]
    omults.append(o[-1] - oinds[-1] + 1)

    n = len(o)
    sn = sum(omults)
    if n != sn:
        raise ValueError("Something Went Wrong With O Sites")
    return osites, omults, oinds


def read_cif(cif):
    atoms, _ = fix_cif(cif)
    ts, tm, tinds = get_tsites(cif)
    os_, om, oinds = get_osites(cif)
    return atoms, ts, tm, tinds, os_, om, oinds


def cif_site_labels(cif):
    _, ts, tm, tinds, os_, om, oinds = read_cif(cif)
    labels = {}
    for i, t in enumerate(ts):
        for j in range(tm[i]):
            labels[tinds[i] + j] = t

    for i, o in enumerate(os_):
        for j in range(om[i]):
            labels[oinds[i] + j] = o

    return labels
