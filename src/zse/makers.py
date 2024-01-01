from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.atoms import Atoms

from zse.cation import monovalent
from zse.collections.framework import get_framework
from zse.substitute import tsub
from zse.t_utilities import get_ratio, get_T_info
from zse.utilities import site_labels

if TYPE_CHECKING:
    from ase.atoms import Atoms


def make_iza_zeolite(code: str) -> Atoms:
    """
    Make an idealized zeolite from an IZA code, populate the atoms.info
    dictionary with the framework name, and add a labels array to the
    atoms object.
    """
    zeolite = get_framework(code)
    labels = site_labels(zeolite, code)
    zeolite.set_array("labels", np.array(list(labels.values())))
    zeolite.info["framework"] = code
    return zeolite


def make_all_exchanged_zeolites(
    code: str,
    heteroatom: str,
    cation: str,
    ignored_T_indices: list[int] | None = None,
) -> list[Atoms]:
    """
    Enumerate all unique T sites and, for each, exchange a single Si atom with a heteroatom.
    Each is charge balanced with the specified cation, and all heteratom-cation configurations
    are returned. Indices in `ignored_T_indices` will not be substituted. If `min_heteroatom_dist`
    is specified, structures with inter-heteroatom distances less than this value will be discarded.

    Limitations:
    - Only supports monovalent cations currently.
    """

    zeolite = make_iza_zeolite(code)
    zeolite.info["heteroatom"] = heteroatom
    zeolite.info["cation"] = cation

    T_info = get_T_info(zeolite, code, ignored_T_indices=ignored_T_indices)

    zeolites = []
    for T_label, T_indices in T_info.items():
        T_index = T_indices[0]
        tsubbed_zeolite = tsub(zeolite, T_index, heteroatom)
        exchanged_zeolites, ring_locations = monovalent(
            tsubbed_zeolite, T_index, cation
        )
        for j, exchanged_zeolite in enumerate(exchanged_zeolites):
            exchanged_zeolite.info["Si_heteroatom_ratio"] = get_ratio(
                exchanged_zeolite, heteroatom=heteroatom
            )
            exchanged_zeolite.info["heteroatom_T_site"] = T_label
            exchanged_zeolite.info["heteroatom_index"] = T_index
            exchanged_zeolite.info["cation_ring"] = ring_locations[j]
            exchanged_zeolite.info["cation_index"] = len(exchanged_zeolite) - 1

            zeolites.append(exchanged_zeolite)
    return zeolites
