__all__ = ["tsub", "nest"]

from copy import deepcopy
from typing import TYPE_CHECKING

from ase import neighborlist
from ase.atoms import Atoms
from ase.build import molecule

from zse.cation import monovalent
from zse.substitute import get_ratio, tsub
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
    Enumerate all unique T sites and, for each, exchange a single Si atom with a heteroatom. Multiple
    configurations for the heteroatom are considered, and all are returned.

    Only supports monovalent cations currently.
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
        T_index = T_indices[
            0
        ]  # TODO: we want to consider all unique T sites; this is only true for all-Si zeolites.
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


def tsub(atoms, index, new_atom):
    """
    atoms should be an ase atoms object
    index is the index of the atom(s) you would like to substitute
    new_atom is the elemental symbol of the atom you want to replace index with.
    """
    z = deepcopy(atoms)
    symbols = z.get_chemical_symbols()
    if isinstance(index, int):
        index = [index]
    for i in index:
        symbols[i] = new_atom

    z.set_chemical_symbols(symbols)

    return z


def nest(atoms, index):
    z = deepcopy(atoms)

    position = z[index].position  # position of that t site

    # This centers the atom object on the T site we want to remove
    center = z.get_center_of_mass()
    trans = center - position
    z.translate(trans)
    z.wrap()

    # get the neighbor list
    cutoff = neighborlist.natural_cutoffs(z, mult=1.05)
    nl = neighborlist.NeighborList(
        cutoffs=cutoff, self_interaction=False, bothways=True
    )
    nl.update(z)
    oxygens = nl.get_neighbors(index)[0]

    # add a hydrogen next to each neighbor oxygen

    for o in oxygens:
        vector = position - z[o].position
        hyd = molecule("H")
        new_location = z[o].position + vector / (vector**2).sum() ** 0.5
        hyd.translate(new_location)
        z = z + hyd

    # recenter the atoms back to their original position and delete the Si
    z.translate(-trans)
    z.wrap()
    del z[index]
    return z
