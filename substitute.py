from ase.io import read, write

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
