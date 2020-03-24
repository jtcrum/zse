from ase.io import read, write
from ase import Atoms, Atoms
import math
import os.path
import numpy as np
import os
from ase.build import molecule

def isolated(atoms,tatom,path = '.'):
    '''
    atoms should be an ase atoms object
    tatom should be the index of the hetero atom with which you want to place
    hydrogen atoms around
    path is the location you want to save the new structures
    new structures are saved with the format of:
    D-INDEX OF OXYGEN THAT H IS BOUND TO
    '''

    lattice = atoms.copy()
    total_oxygen = [atom.index for atom in atoms if atom.symbol == 'O']
    tsites = [tatom]
    total_silicon = [atom.index for atom in atoms if atom.symbol == 'Si']

    n = len(tsites)

    hydrogen = []
    for i in range(n):
        hydrogen.append(len(atoms)+i)

    oxygens = []
    for j in tsites:
        tmp = []
        for k in total_oxygen:
            distance = lattice.get_distance(j,k,mic=True)
            if distance < 2.0:
                tmp.append(k)
                if len(tmp) == 4:
                    oxygens.append(tmp)
                    tmp = []
    oxygens = np.array(oxygens)
    first_oxygens = oxygens[0]

    silicons = []
    for l in first_oxygens:
        tmp = []
        for m in total_silicon:
            distance = lattice.get_distance(l,m,mic=True)
            if distance < 2.0:
                silicons.append(m)
    silicons = np.array(silicons)

    adsorbate = molecule('H')
    adsorbate.translate([0,0,0])
    H_lattice = lattice + adsorbate

    cwd = os.getcwd()
    traj = []
    for l in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[tsites]
        H_lattice.translate(diff)
        H_lattice.wrap()
        H_lattice.set_distance(first_oxygens[l], hydrogen[0], 0.98, fix = 0)
        H_lattice.set_angle(int(tsites[0]), first_oxygens[l], hydrogen[0], 109.6, mask = None)
        H_lattice.set_angle(int(silicons[l]), first_oxygens[l],hydrogen[0], 109.6, mask = None)
        H_lattice.set_dihedral(int(tsites[0]), first_oxygens[l], silicons[l], hydrogen[0], 180, mask = None)
        H_lattice.translate(-1*diff)
        H_lattice.wrap()
        traj += [Atoms(H_lattice)]

        if not os.path.exists(path+'/D-' + str(first_oxygens[l])):
            os.makedirs(path+'/D-' + str(first_oxygens[l]))
        write(path+'/D-'  + str(first_oxygens[l]) +'/POSCAR',H_lattice, sort = True)

    return traj , oxygens

def paired(atoms, tatoms, path = '.'):
    '''
    atoms should be an ase atoms object
    tatoms should be the indices of the hetero atoms with which you want to place
    hydrogen atoms around
    alternatively, tatoms can be a symbol i.e. 'Al' if both hetero atoms are the same element
    path is the location you want to save the new structures
    new structures are saved with the format of:
    D-INDEX OF OXYGEN THAT H IS BOUND TO
    '''
    lattice = atoms.copy()
    total_oxygen = [atom.index for atom in atoms if atom.symbol == 'O']

    if isinstance(tatoms,list):
        tsites = tatoms
    elif isinstance(tatoms, str):
        tsites = [atom.index for atom in atoms if atom.symbol == tatoms]

    total_silicon = [atom.index for atom in atoms if atom.symbol == 'Si']
    n = len(atoms)
    Hydrogen = []
    for i in range(n):
        Hydrogen.append(len(atoms)+i)

    oxygens = []
    for j in tsites:
        tmp = []
        for k in total_oxygen:
            distance = lattice.get_distance(j,k,mic=True)
            if distance < 2.0:
                tmp.append(k)
                if len(tmp) == 4:
                    oxygens.append(tmp)
                    tmp = []
    oxygens = np.array(oxygens)
    first_oxygens = oxygens[0]
    second_oxygens = oxygens[1]

    silicons = []
    for l in oxygens[0]:
        tmp = []
        for m in total_silicon:
            distance = lattice.get_distance(l,m,mic=True)
            if distance < 2.0:
                silicons.append(m)
    silicons = np.array(silicons)
    first_silicons = silicons

    silicons = []
    for l in oxygens[1]:
        tmp = []
        for m in total_silicon:
            distance = lattice.get_distance(l,m,mic=True)
            if distance < 2.0:
                silicons.append(m)
    silicons = np.array(silicons)
    second_silicons = silicons

    adsorbate = molecule('H')
    adsorbate.translate([0, 0, 0])
    H_lattice = lattice + adsorbate + adsorbate

    cwd = os.getcwd()

    first_oxygens=first_oxygens.astype(int)
    first_silicons = first_silicons.astype(int)
    tsites = np.array(tsites)
    tsites = tsites.astype(int)
    Hydrogen = np.array(Hydrogen)
    Hydrogen = Hydrogen.astype(int)
    traj = []
    for l in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[tsites[0]]
        H_lattice.translate(diff)
        H_lattice.wrap()
        first_distance = H_lattice.set_distance(first_oxygens[l], Hydrogen[0], 0.98, fix=0)
        first_Al_angle = H_lattice.set_angle(int(tsites[0]), first_oxygens[l], Hydrogen[0], 109.6, mask=None)
        first_Si_angle = H_lattice.set_angle(int(first_silicons[l]), first_oxygens[l], Hydrogen[0], 109.6, mask=None)
        first_dihedral = H_lattice.set_dihedral(int(tsites[0]), first_oxygens[l], first_silicons[l], Hydrogen[0], 180, mask=None)
        H_lattice.translate(-1*diff)
        H_lattice.wrap()
        for k in range(4):
            center = H_lattice.get_center_of_mass()
            positions = atoms.get_positions()
            diff = center - positions[tsites[1]]
            H_lattice.translate(diff)
            H_lattice.wrap()
            second_distance = H_lattice.set_distance(second_oxygens[k], Hydrogen[1], 0.98, fix=0)
            second_Al_angle = H_lattice.set_angle(int(tsites[1]), second_oxygens[k], Hydrogen[1], 109.6, mask=None)
            second_Si_angle = H_lattice.set_angle(int(second_silicons[k]), second_oxygens[k], Hydrogen[1], 109.6, mask=None)
            second_dihedral = H_lattice.set_dihedral(int(tsites[1]), second_oxygens[k], second_silicons[k], Hydrogen[1], 180, mask=None)
            H_lattice.translate(-1*diff)
            H_lattice.wrap()
            # makes directory of each First_Oxygen-Second_Oxygen case with terminal hydrogen
            if not os.path.exists('{0}/D-{1}-{2}'.format(path,str(first_oxygens[l]),str(second_oxygens[k]))):
                os.makedirs('{0}/D-{1}-{2}'.format(path,str(first_oxygens[l]),str(second_oxygens[k])))

            write('{0}/D-{1}-{2}/POSCAR'.format(path,str(first_oxygens[l]),str(second_oxygens[k])),H_lattice, sort = True)
            traj += [Atoms(H_lattice)]

    return traj , oxygens
