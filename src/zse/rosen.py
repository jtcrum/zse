"""Scattered utilities for Andrew Rosen's zeolite projects."""
from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms

from zse.cation import monovalent
from zse.collections import framework
from zse.substitute import tsub
from zse.utilities import site_labels

if TYPE_CHECKING:
    from typing import Any


def make_iza_zeolite(code: str) -> Atoms:
    """
    Make an idealized zeolite from an IZA code, populate the atoms.info
    dictionary with the framework name, and add a labels array to the
    atoms object.
    """
    zeolite = framework(code)
    labels = site_labels(zeolite, code)
    zeolite.set_array("labels", np.array(list(labels.values())))
    zeolite.info["framework"] = code
    return zeolite

def get_ratio(atoms: Atoms, heteroatom: str, offset: int = 0) -> float:
    """
    Calculate the Si/heteroatom ratio of a zeolite.
    """
    n_Si = len([atom for atom in atoms if atom.symbol == "Si"]) - offset
    n_heteroatom = len([atom for atom in atoms if atom.symbol == heteroatom]) + offset

    return n_Si / n_heteroatom


def _prep_labels(d: dict, labels: list[Any]) -> None:
    """
    Prepare the atoms.info entries
    """
    for label in labels:
        if label not in d:
            d[label] = []


def get_T_info(
    zeolite: Atoms, code: str, ignored_T_indices: list[int] = None
) -> dict[str, list[int]]:
    """
    Get T-site info for a zeolite. Returns a dictionary of the form
    {"T1": [0, 1], "T2": [1, 2], ...} where the keys are the T-site labels and the values
    are the indices of the T-site in the zeolite. If ignored_indices is
    specified, then these will be excluded from the returned indices.
    """
    ignored_T_indices = ignored_T_indices or []

    labels = list(site_labels(zeolite, code).values())
    unique_T_labels = np.unique([T for T in labels if "T" in T]).tolist()
    T_info = {}

    for T_label in unique_T_labels:
        T_indices = [
            i
            for i, label in enumerate(labels)
            if label == T_label and i not in ignored_T_indices
        ]
        T_info[T_label] = T_indices
    return T_info


def get_min_T_distance(atoms: Atoms, T_symbols: str | list[str]) -> float:
    """
    Get the minimum distance between all heteroatom pairs in a zeolite.
    """
    if isinstance(T_symbols, str):
        T_symbols = [T_symbols]
    heteroatom_sites = [atom.index for atom in atoms if atom.symbol in T_symbols]
    if len(heteroatom_sites) > 1:
        heteroatom_positions = atoms[heteroatom_sites].get_all_distances(mic=True)
        return np.min(heteroatom_positions[heteroatom_positions > 0])
    else:
        return np.inf

def exchange_unique_T_sites(
    zeolite: Atoms,
    code: str,
    heteroatom: str,
    cation: str,
    ignored_T_indices: list[int] | None = None,
    min_heteroatom_dist: float | None = 3.5,
) -> list[Atoms]:
    """
    Enumerate all unique T sites and, for each, exchange a single Si atom with a heteroatom. Multiple
    configurations for the heteroatom are considered, and all are returned.

    Only supports monovalent cations currently.
    """
    zeolite = deepcopy(zeolite)
    zeolites = []

    zeolite.info["framework"] = code
    zeolite.info["heteroatom"] = heteroatom
    zeolite.info["cation"] = cation

    T_info = get_T_info(zeolite, code, ignored_T_indices=ignored_T_indices)

    for T_label, T_indices in T_info.items():
        T_index = T_indices[0] # TODO: we want to consider all unique T sites; this is only true for all-Si zeolites.
        exchanged_zeolites, ring_locations = monovalent(
            tsub(zeolite, T_index, heteroatom), T_index, cation
        )
        for j, exchanged_zeolite_ in enumerate(exchanged_zeolites):
            exchanged_zeolite = deepcopy(exchanged_zeolite_)

            _prep_labels(
                exchanged_zeolite.info,
                [
                    "heteroatom_T_sites",
                    "heteroatom_indices",
                    "cation_rings",
                    "cation_indices",
                    "cation_CrystalNN",
                ],
            )

            exchanged_zeolite.info["Si_heteroatom_ratio"] = get_ratio(
                exchanged_zeolite, heteroatom
            )
            exchanged_zeolite.info["heteroatom_T_sites"].append(T_label)
            exchanged_zeolite.info["heteroatom_indices"].append(T_index)
            exchanged_zeolite.info["cation_rings"].append(ring_locations[j])
            exchanged_zeolite.info["cation_indices"].append(len(exchanged_zeolite) - 1)
            exchanged_zeolite.info["cation"] = cation

            if min_heteroatom_dist:
                min_dist = get_min_T_distance(exchanged_zeolite, heteroatom)
                if min_dist < min_heteroatom_dist:
                    continue
            zeolites.append(exchanged_zeolite)
    return zeolites
