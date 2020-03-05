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
    tatom should be the symbol of the hetero atom with which you want to place
    hydrogen atoms around
    path is the location you want to save the new structures
    new structures are saved with the format of:
    D-INDEX OF OXYGEN THAT H IS BOUND TO
    '''

    lattice = atoms.copy()
    total_oxygen = [atom.index for atom in atoms if atom.symbol == 'O']
    tsites = [atom.index for atom in atoms if atom.symbol == tatom]
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
