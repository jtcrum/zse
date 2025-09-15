__all__ = ["tsub", "nest"]

from ase.io import read, write
import numpy as np
from ase import neighborlist
from ase.build import molecule
from itertools import permutations
from ase import Atoms


def tsub(atoms, index, new_atom):
    """
    atoms should be an ase atoms object
    index is the index of the atom(s) you would like to substitute
    new_atom is the elemental symbol of the atom you want to replace index with.
    """
    z = atoms.copy()
    symbols = z.get_chemical_symbols()
    if isinstance(index, int):
        index = [index]
    for i in index:
        symbols[i] = new_atom

    z.set_chemical_symbols(symbols)

    return z


def nest(z, index):
    """
    z should be an ase atoms object (must be a t-site)
    index is the index of the atom you would like to delete and replace with a defect
    """
    position = z[index].position  # position of that t site

    # This centers the atom object on the T site we want to remove
    center = z.get_center_of_mass()
    trans = center - position
    z.translate(trans)
    z.wrap()

    # get the neighbor list
    cutoff = neighborlist.natural_cutoffs(z, mult=1.05)
    nl = neighborlist.NeighborList(cutoffs=cutoff, self_interaction=False, bothways=True)
    nl.update(z)
    oxygens = nl.get_neighbors(index)[0]

    if len(oxygens) == 2:
        try:
            raise SystemExit("error in code want to exit")
        except:
            print("program is still open")

    # iterate over proton configurations
    o_list = list(permutations(oxygens, 4))[
        0:6
    ]  # permutations are cyclic so we only need the first 6 values

    # for each cycle, add a proton poitning towards the next oxygen in the cycle
    nest_traj = []  # store structures
    for cycle in o_list:
        # O ids
        atoms = z.copy()
        A = cycle[0]
        B = cycle[1]
        C = cycle[2]
        D = cycle[3]
        # O positions
        pA = atoms[A].position
        pB = atoms[B].position
        pC = atoms[C].position
        pD = atoms[D].position
        # Vectors between O to place H
        vA = (
            atoms.get_distance(A, B, vector=True, mic=True) / 2.5
        )  # vetor is divided by 2.5 so its not eactly between two protons
        vB = (
            atoms.get_distance(B, C, vector=True, mic=True) / 2.5
        )  # closer to the first O atom in each vector
        vC = atoms.get_distance(C, D, vector=True, mic=True) / 2.5
        vD = atoms.get_distance(D, A, vector=True, mic=True) / 2.5

        # Add H at each position in cycle
        H_atoms = Atoms("H4", positions=[pA + vA, pB + vB, pC + vC, pD + vD])

        new_atoms = atoms + H_atoms

        new_atoms.translate(-trans)
        new_atoms.wrap()
        del new_atoms[index]

        nest_traj.append(new_atoms)

    return nest_traj
