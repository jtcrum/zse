"""Tools for determining uniqueness"""
from __future__ import annotations
from scipy.spatial.distance import pdist, squareform

import numpy as np
from ase.atoms import Atoms
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pymatgen.core.structure import Structure
def get_soap_distances(
    atoms: Atoms, indices: list[int], rcut: float = 6.0, nmax: int = 8, lmax: int = 6
) -> np.ndarray:
    """
    Get the SOAP distance between all specified indices in a zeolite.
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

def get_unique_structures(zeolites: list[Atoms]) -> list[Structure] | list[Atoms]:
    """
    Get a unique list of structures from a list of Atoms objects.
    """

    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.io.ase import AseAtomsAdaptor
    structures = [AseAtomsAdaptor().get_structure(atoms) for atoms in zeolites]
    unique_structures = [s[0] for s in StructureMatcher().group_structures(structures)]
    return unique_structures
