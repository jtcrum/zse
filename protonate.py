from ase.io import read, write
from ase import Atoms, Atoms
import math
import os.path
import numpy as np
import os
import shutil
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
            os.mkdir(path+'/D-' + str(first_oxygens[l]))
        write(path+'/D-'  + str(first_oxygens[l]) +'/POSCAR',H_lattice, sort = True)

    return traj

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
    total_oxygen = [atom.index for atom in atoms if atom.symnbol == 'O']

    if isinstance(tatoms,list):
        tsites = tatoms
    elif isinstance(tatoms, str):
        tsites = [atom.index for atom in atoms if atom.symmbol == tatoms]

    total_silicon = [atom.index for atom in atoms if atom.symbol == 'Si']
    n = len(atoms)
    Hydrogen = []
    for i in range(n):
        Hydorgen.append(len(atoms)+i)

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

    First_Oxygens=First_Oxygens.astype(int)
    First_Silicons = First_Silicons.astype(int)
    Aluminium = np.array(Aluminium)
    Aluminium = Aluminium.astype(int)
    Hydrogen = np.array(Hydrogen)
    Hydrogen = Hydrogen.astype(int)
    traj = []
    for l in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[Aluminium[0]]
        H_lattice.translate(diff)
        H_lattice.wrap()
        first_distance = H_lattice.set_distance(First_Oxygens[l], Hydrogen[0], 0.98, fix=0)
        first_Al_angle = H_lattice.set_angle(int(Aluminium[0]), First_Oxygens[l], Hydrogen[0], 109.6, mask=None)
        first_Si_angle = H_lattice.set_angle(int(First_Silicons[l]), First_Oxygens[l], Hydrogen[0], 109.6, mask=None)
        first_dihedral = H_lattice.set_dihedral(int(Aluminium[0]), First_Oxygens[l], First_Silicons[l], Hydrogen[0], 180, mask=None)
        H_lattice.translate(-1*diff)
        H_lattice.wrap()
        for k in range(4):
            center = H_lattice.get_center_of_mass()
            positions = atoms.get_positions()
            diff = center - positions[Aluminium[1]]
            H_lattice.translate(diff)
            H_lattice.wrap()
            second_distance = H_lattice.set_distance(Second_Oxygens[k], Hydrogen[1], 0.98, fix=0)
            second_Al_angle = H_lattice.set_angle(int(Aluminium[1]), Second_Oxygens[k], Hydrogen[1], 109.6, mask=None)
            second_Si_angle = H_lattice.set_angle(int(Second_Silicons[k]), Second_Oxygens[k], Hydrogen[1], 109.6, mask=None)
            second_dihedral = H_lattice.set_dihedral(int(Aluminium[1]), Second_Oxygens[k], Second_Silicons[k], Hydrogen[1], 180, mask=None)
            H_lattice.translate(-1*diff)
            H_lattice.wrap()
            # makes directory of each First_Oxygen-Second_Oxygen case with terminal hydrogen
            if not os.path.exists('{0}/D-{1}-{2}'.format(path,str(first_oxygens[l]),str(second_oxygens[k]))):
                os.mkdir('{0}/D-{1}-{2}'.format(path,str(first_oxygens[l]),str(second_oxygens[k])))

            write('{0}/D-{1}-{2}/POSCAR'.format(path,str(first_oxygens[l]),str(second_oxygens[k])),H_lattice, sort = True)
            traj += [Atoms(H_lattice)]
