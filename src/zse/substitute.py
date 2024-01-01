from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING

from ase import neighborlist
from ase.atoms import Atoms
from ase.build import molecule

from zse.cation import monovalent
from zse.t_utilities import get_ratio, get_T_info
from zse.utilities import make_iza_zeolite

if TYPE_CHECKING:
    from typing import Any


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


def exchange_unique_T_sites(
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
            exchanged_zeolite.info["cation_index"] = (len(exchanged_zeolite) - 1)
            exchanged_zeolite.info["cation"] = cation

            zeolites.append(exchanged_zeolite)
    return zeolites
