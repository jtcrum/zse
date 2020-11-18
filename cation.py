from ase.io import read, write
from ase import Atoms, Atoms
import math
import os.path
import numpy as np
import os
from ase.build import molecule

def divalent(atoms,M,path = None):
    '''
    This function will place one divalent cation into a zeolite framework that
    contains two negative charge centers.
    The cation is placed at each of 6 rings around each of the aluminum atoms.
    This creates 12 structures. Depending on the placement of your Al, some of
    the strucutres may not be unique.

    The structures will be placed in the path provided as an input.
    Format of the structures folder names are:
    D-aluminum_index-oxygen1_index,oxygen2_index.

    The original version of this code was written by Sichi Li in 2017.
    '''

    total_oxygen = [atom.index for atom in atoms if atom.symbol == 'O']
    aluminum = [atom.index for atom in atoms if atom.symbol == 'Al']
    cation = [len(atoms)]

    oxygens = []
    for j in aluminum:
        tmp = []
        for k in total_oxygen:
            distance = atoms.get_distance(k,j,mic=True)
            if distance < 1.7:
                tmp.append(k)
                if len(tmp)  == 4:
                    oxygens.append(tmp)
                    tmp = []
    oxygens = np.array(oxygens)
    traj = []
    for l in range(len(aluminum)):
        for i in range(len(oxygens[l,:])):
            for j in (x for x in range(len(oxygens[l,:])) if x > i):
                totpos = atoms.get_positions()
                a = totpos[oxygens[l,i],:]
                b = totpos[oxygens[l,j],:]
                c = totpos[aluminum[l],:]
                pos = a-c+b
                tpos = pos.reshape(1,3)
                adsorbate = Atoms(M)
                adsorbate.set_positions(tpos)
                M_lattice = atoms+adsorbate
                traj+=[M_lattice]

                if path:
                    os.makedirs('{0}/D-{1}-{2}-{3}'.format(path,str(aluminum[l]),str(oxygens[l,j]),str(oxygens[l,i])))

                    write('{0}/D-{1}-{2}-{3}/POSCAR'.format(path,str(aluminum[l]),str(oxygens[l,j]),str(oxygens[l,i])),M_lattice, sort = True)
    return traj

def monovalent_old(atoms,symbol,path=None):

    '''
    This function will place one monovalent cation into a zeolite framework that
    contains one aluminum atom.
    The cation is placed at each of 6 rings around the aluminum atoms.
    This creates 6 structures.

    The structures will be placed in the path provided as an input.
    Format of the structures folder names are:
    D-aluminum_index-oxygen1_index,oxygen2_index.
    '''

    total_oxygen = [atom.index for atom in atoms if atom.symbol == 'O']
    aluminum = [atom.index for atom in atoms if atom.symbol == 'Al']
    cation = [len(atoms)]
    oxygens = []
    for j in aluminum:
        tmp = []
        for k in total_oxygen:
            distance = atoms.get_distance(k,j,mic=True)
            if distance < 2:
                tmp.append(k)
                if len(tmp)  == 4:
                    oxygens.append(tmp)
                    tmp = []

    oxygens = np.array(oxygens)
    traj = []

    for l in range(len(aluminum)):
        for i in range(len(oxygens[l,:])):
            for j in (x for x in range(len(oxygens[l,:])) if x > i):
                totpos = atoms.get_positions()
                a = totpos[oxygens[l,i],:]
                b = totpos[oxygens[l,j],:]
                c = totpos[aluminum[l],:]
                pos = a-c+b
                tpos = pos.reshape(1,3)
                adsorbate = Atoms(symbol)
                adsorbate.set_positions(tpos)
                M_lattice = atoms+adsorbate
                traj+=[M_lattice]

                if path:
                    os.makedirs('{0}/D-{1}-{2}-{3}'.format(path,str(aluminum[l]),str(oxygens[l,j]),str(oxygens[l,i])),exist_ok=True)

                    write('{0}/D-{1}-{2}-{3}/POSCAR'.format(path,str(aluminum[l]),str(oxygens[l,j]),str(oxygens[l,i])),M_lattice, sort = True)

    return traj

def monovalent(atoms,index,symbol,framework,included_rings=None,path=None):

    from zse import rings
    import networkx as nx
    from ase import neighborlist
    from ase.data import covalent_radii, chemical_symbols

    radii = {chemical_symbols[i]: covalent_radii[i] for i in range(len(chemical_symbols))}

    possible = rings.get_fwrings(framework)
    possible = possible*2
    if included_rings == None:
        included_rings = []
        for p in possible:
            if p > 8:
                included_rings.append(p)
    atoms2 = atoms.copy()
    cell = atoms2.get_cell_lengths_and_angles()[:3]
    repeat = []
    maxring = max(possible)
    for i,c in enumerate(cell):
        if c/2 < maxring/2+5:
            l = c
            re = 1
            while l/2 < maxring/2+5:
                re +=1
                l = c*re

            repeat.append(re)
        else:
            repeat.append(1)
    atoms2 = atoms2.repeat(repeat)
    center = atoms2.get_center_of_mass()
    trans = center - atoms2.positions[index]
    atoms2.translate(trans)
    atoms2.wrap()

    cutoff = neighborlist.natural_cutoffs(atoms2, mult = 0.95)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms2)
    matrix = nl.get_connectivity_matrix(sparse = False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)
    r = rings.find_o_rings(G,index,possible)
    paths = rings.remove_dups(r)
    paths = rings.remove_sec(paths)

    Class = []
    for p in paths:
        Class.append(int(len(p)/2))

    paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
    Class.sort(reverse=True)

    if path:
        class_count = []
        for i in range(len(Class)):
            if i == 0:
                class_count.append(1)
            else:
                counter = 1
                for j in range(i-1):
                    if Class[i]==Class[j]:
                        counter+=1
                class_count.append(counter)


    # get center of mass for each ring, place ion in that ring
    traj= []
    for i in range(len(paths)):
        p = paths[i]
        if len(p) in included_rings:
            positions = atoms2[p].positions
            co = sum(positions)/len(positions)
            co -=trans
            vector = co-atoms[index].position
            vector = np.array(vector)
            vhat = vector/np.linalg.norm(vector)

            bond_length = radii[symbol]+2/radii[symbol]

            new_p = [atoms[index].position +bond_length*vhat]

            adsorbate = Atoms(symbol)
            adsorbate.set_positions(new_p)
            c_atoms = atoms+adsorbate
            c_atoms.wrap()
            traj+=[c_atoms]

            # write POSCAR for each structure

            if path:
                os.makedirs('{0}/D-{1}MR-{2}'.format(path,str(int(len(p)/2)),str(class_count[i])),exist_ok=True)

                write('{0}/D-{1}MR-{2}/POSCAR'.format(path,str(int(len(p)/2)),str(class_count[i])),c_atoms, sort = True)


    return Class ,paths, traj, repeat
