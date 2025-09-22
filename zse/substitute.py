"""Tools to substitute atoms in an ASE atoms object
to create defects or doped structures.
"""

from itertools import permutations

from ase import Atoms, neighborlist

__all__ = ["nest", "tsub"]


def tsub(atoms: Atoms, index: int | list[int], new_atom: str) -> Atoms:
    """Substitute atom(s) in an ASE atoms object with a new element.

    Args:
        atoms (Atoms): The ASE atoms object to modify.
        index (int | list[int]): The index or list of indices of the atom(s) to substitute.
        new_atom (str): The elemental symbol of the new atom to substitute in.

    Returns:
        Atoms: A new ASE atoms object with the specified atom(s) substituted.
    """
    z = atoms.copy()
    symbols = z.get_chemical_symbols()
    if isinstance(index, int):
        index = [index]
    for i in index:
        symbols[i] = new_atom

    z.set_chemical_symbols(symbols)

    return z


def nest(z: Atoms, index: int) -> list[Atoms]:
    """Create a list of structures with a nest defect at a specified T-site.

    Args:
        z (Atoms): The ASE atoms object representing the zeolite framework.
        index (int): The index of the T-site atom to remove and create the nest defect
    Returns:
        list[Atoms]: A list of ASE atoms objects, each representing a different
            configuration of the nest defect.
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
        except Exception as e:
            print(f"Program is still open!! Error message: {e}")

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
