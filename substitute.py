__all__ = ['tsub','nest']

from ase.io import read, write
import numpy as np
from ase import neighborlist
from ase.build import molecule


def tsub(atoms,index,new_atom):
    '''
    atoms should be an ase atoms object
    index is the index of the atom(s) you would like to substitute
    new_atom is the elemental symbol of the atom you want to replace index with.
    '''

    symbols = atoms.get_chemical_symbols()
    if isinstance(index, int):
        index = [index]
    for i in index:
        symbols[i]=new_atom

    atoms.set_chemical_symbols(symbols)

    return atoms

def nest(atoms,index):

    position = atoms[index].position # position of that t site

    # This centers the atom object on the T site we want to remove
    center = atoms.get_center_of_mass()
    trans = center-position
    atoms.translate(trans)
    atoms.wrap()

    # get the neighbor list
    cutoff = neighborlist.natural_cutoffs(atoms,mult=1.05)
    nl = neighborlist.NeighborList(cutoffs=cutoff, self_interaction = False, bothways = True)
    nl.update(atoms)
    oxygens = nl.get_neighbors(index)[0]

    # add a hydrogen next to each neighbor oxygen

    for o in oxygens:
        vector = position - atoms[o].position
        hyd = molecule('H')
        new_location = atoms[o].position + vector / (vector**2).sum()**.5
        hyd.translate(new_location)
        atoms = atoms + hyd


    # recenter the atoms back to their original position and delete the Si
    atoms.translate(-trans)
    atoms.wrap()
    del(atoms[index])
    return atoms
