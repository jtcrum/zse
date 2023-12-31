"""Get info about T sites in a zeolite."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from zse.utilities import site_labels

if TYPE_CHECKING:
    from ase.atoms import Atoms


def get_ratio(atoms: Atoms, heteroatom: str) -> float:
    """
    Calculate the Si/heteroatom ratio of a zeolite.
    """
    n_Si = len([atom for atom in atoms if atom.symbol == "Si"])
    n_heteroatom = len([atom for atom in atoms if atom.symbol == heteroatom])

    return n_Si / n_heteroatom


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
