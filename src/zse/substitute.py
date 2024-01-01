from __future__ import annotations

from copy import deepcopy

from ase import neighborlist
from ase.build import molecule


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
