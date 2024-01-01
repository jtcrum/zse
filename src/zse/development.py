from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

from ase.atoms import Atoms

from zse.cation import monovalent
from zse.substitute import tsub
from zse.t_utilities import get_min_T_distance, get_ratio, get_T_info

if TYPE_CHECKING:
    from typing import Any


def exchange_unique_T_sites(
    zeolite: Atoms,
    code: str,
    heteroatom: str,
    cation: str,
    ignored_T_indices: list[int] | None = None,
    min_heteroatom_dist: float | None = 3.5,
) -> list[Atoms]:
    """
    Enumerate all unique T sites and, for each, exchange a single Si atom with a heteroatom.
    Each is charge balanced with the specified cation, and all heteratom-cation configurations
    are returned. Indices in `ignored_T_indices` will not be substituted. If `min_heteroatom_dist`
    is specified, structures with inter-heteroatom distances less than this value will be discarded.

    Limitations:
    - Unique T sites are determined without consideration of the underlying atom identity.
    - Only supports monovalent cations currently.
    """

    def _prep_labels(d: dict, labels: list[Any]) -> None:
        """
        Prepare the atoms.info entries
        """
        for label in labels:
            if label not in d:
                d[label] = []

    zeolite = deepcopy(zeolite)
    zeolites = []

    zeolite.info["framework"] = code
    zeolite.info["heteroatom"] = heteroatom
    zeolite.info["cation"] = cation

    T_info = get_T_info(zeolite, code, ignored_T_indices=ignored_T_indices)

    for T_label, T_indices in T_info.items():
        T_index = T_indices[0]
        tsubbed_zeolite = tsub(zeolite, T_index, heteroatom)
        exchanged_zeolites, ring_locations = monovalent(
            tsubbed_zeolite, T_index, cation
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
                exchanged_zeolite, heteroatom=heteroatom
            )
            exchanged_zeolite.info["heteroatom_T_sites"].append(T_label)
            exchanged_zeolite.info["heteroatom_indices"].append(T_index)
            exchanged_zeolite.info["cation_rings"].append(ring_locations[j])
            exchanged_zeolite.info["cation_indices"].append(len(exchanged_zeolite) - 1)
            exchanged_zeolite.info["cation"] = cation

            if min_heteroatom_dist:
                min_dist = get_min_T_distance(exchanged_zeolite, T_symbols=heteroatom)
                if min_dist < min_heteroatom_dist:
                    continue
            zeolites.append(exchanged_zeolite)
    return zeolites
