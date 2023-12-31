"""Scattered utilities for Andrew Rosen's zeolite projects."""
from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.spatial.distance import pdist, squareform

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


def get_unique_structures(zeolites: list[Atoms]) -> list[Structure] | list[Atoms]:
    """
    Get a unique list of structures from a list of Atoms objects.
    """
    structures = [AseAtomsAdaptor().get_structure(atoms) for atoms in zeolites]
    unique_structures = [s[0] for s in StructureMatcher().group_structures(structures)]
    return unique_structures


def get_ratio(atoms: Atoms, heteroatom: str, offset: int = 0) -> float:
    """
    Calculate the Si/heteroatom ratio of a zeolite.
    """
    n_Si = len([atom for atom in atoms if atom.symbol == "Si"]) - offset
    n_heteroatom = len([atom for atom in atoms if atom.symbol == heteroatom]) + offset

    return n_Si / n_heteroatom


def prep_labels(d: dict, labels: list[Any]) -> None:
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


def get_min_heteroatom_distance(atoms: Atoms, heteroatom: str) -> float:
    """
    Get the minimum distance between all heteroatom pairs in a zeolite.
    """
    min_dist = np.inf
    heteroatom_sites = [atom.index for atom in atoms if atom.symbol == heteroatom]
    if len(heteroatom_sites) > 1:
        heteroatom_positions = atoms[heteroatom_sites].get_positions()
        for i in range(len(heteroatom_positions)):
            for j in range(i + 1, len(heteroatom_positions)):
                distance = atoms.get_distance(
                    heteroatom_sites[i], heteroatom_sites[j], mic=True
                )
                if distance < min_dist:
                    min_dist = distance
    return min_dist


def get_soap_distances(
    atoms: Atoms, indices: list[int], rcut: float = 6.0, nmax: int = 8, lmax: int = 6
) -> np.ndarray:
    """
    Get the SOAP distance between all heteroatom pairs in a zeolite.
    """
    from dscribe.descriptors import SOAP

    soap = SOAP(
        species=list(set(atoms.get_chemical_symbols())),
        periodic=True,
        rcut=rcut,
        nmax=nmax,
        lmax=lmax,
    )
    soap_values = soap.create(atoms, centers=indices)
    pairwise_distances = pdist(soap_values, "euclidean")
    return squareform(pairwise_distances)


def find_unique_samples(X: np.ndarray, tol: float = 1e-4):
    """
    Find unique samples in a square distance matrix.

    TODO.
    """


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

            prep_labels(
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
                min_dist = get_min_heteroatom_distance(exchanged_zeolite, heteroatom)
                if min_dist < min_heteroatom_dist:
                    continue
            zeolites.append(exchanged_zeolite)
    return zeolites
