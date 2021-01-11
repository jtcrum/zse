__all__ = ['get_os_and_ts','add_one_proton','add_two_protons']

from zse.utilities import *
from ase.io import read,write
from ase import Atoms
from ase.build import molecule
import numpy as np
import os

def get_os_and_ts(atoms,index):
    lattice = atoms.copy()
    total_oxygen = [atom.index for atom in lattice if atom.symbol =='O']
    total_silicon = [atom.index for atom in lattice if atom.symbol =='Si']

    oxygens = []
    for k in total_oxygen:
        distance = lattice.get_distances(index,k,mic=True)
        if distance < 2.0:
            oxygens.append(k)
    oxygens = np.array(oxygens)
    silicons = []
    for l in oxygens:
        tmp = []
        for m in total_silicon:
            distance = lattice.get_distance(l,m,mic=True)
            if distance < 2.0:
                silicons.append(m)
    silicons = np.array(silicons)

    return oxygens, silicons

def add_one_proton(atoms, index, oxygens, silicons, code, path=None):

    labels = site_labels(atoms,code)

    hydrogen = [len(atoms)]

    adsorbate = molecule('H')
    adsorbate.translate([0,0,0])
    H_lattice = atoms + adsorbate

    cwd = os.getcwd()
    traj = []
    locations = []
    for l in range(4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[index]
        H_lattice.translate(diff)
        H_lattice.wrap()
        H_lattice.set_distance(oxygens[l],hydrogen[0],0.98, fix = 0)
        H_lattice.set_angle(int(index),int(oxygens[l]),int(hydrogen[0]), 109.6, mask = None)
        H_lattice.set_angle(int(silicons[l]),int(oxygens[l]),int(hydrogen[0]), 109.6, mask = None)
        H_lattice.set_dihedral(int(index), int(oxygens[l]), int(silicons[l]), hydrogen[0], 180, mask = None)
        H_lattice.translate(-1*diff)
        H_lattice.wrap()
        traj += [Atoms(H_lattice)]
        locations.append(labels[oxygens[l]])

        if path:
            os.makedirs('{0}/D-{1}'.format(path,labels[oxygens[l]]),exist_ok=True)

            write('{0}/D-{1}/POSCAR'.format(path,labels[oxygens[l]]),H_lattice, sort = False)

    return traj, locations

def add_two_protons(atoms, indices, oxygens, silicons, code, path = None):

    labels = site_labels(atoms,code)

    hydrogen =[len(atoms),len(atoms)+1]
    adsorbate = molecule('H')
    H_lattice = atoms + adsorbate + adsorbate

    cwd = os.getcwd()
    locations=[]
    traj = []
    for l in range (4):
        center = H_lattice.get_center_of_mass()
        positions = atoms.get_positions()
        diff = center - positions[indices[0]]
        H_lattice.translate(diff)
        H_lattice.wrap()
        first_distance = H_lattice.set_distance(oxygens[0][l],hydrogen[0], 0.98, fix=0)
        first_Al_angle = H_lattice.set_angle(int(indices[0]), oxygens[0][l], hydrogen[0], 109.6, mask=None)
        first_Si_angle = H_lattice.set_angle(int(silicons[0][l]),oxygens[0][l], hydrogen[0], 109.6, mask = None)
        first_dihedral = H_lattice.set_dihedral(int(indices[0]),oxygens[0][l], silicons[0][l], hydrogen[0], 180, mask = None)
        H_lattice.translate(-1*diff)
        H_lattice.wrap()
        for k in range(4):
            center = H_lattice.get_center_of_mass()
            positions = atoms.get_positions()
            diff = center - positions[indices[1]]
            H_lattice.translate(diff)
            H_lattice.wrap()
            second_distance = H_lattice.set_distance(oxygens[1][k],hydrogen[1],0.98,fix=0)
            second_Al_angle = H_lattice.set_angle(int(indices[1]),oxygens[1][k],hydrogen[1],109.6,mask=None)
            second_Si_angle = H_lattice.set_angle(int(silicons[1][k]),oxygens[1][k],hydrogen[1],109.6, mask=None)
            second_dihedral = H_lattice.set_dihedral(int(indices[1]),oxygens[1][k],silicons[1][k],hydrogen[1],180,mask=None)
            H_lattice.translate(-1*diff)
            H_lattice.wrap()

            traj += [Atoms(H_lattice)]
            locations.append('{0}-{1}'.format(labels[oxygens[0][l]],labels[oxygens[1][k]]))
            if path:
                os.makedirs('{0}/D-{1}-{2}'.format(path,labels[oxygens[0][l]],labels[oxygens[1][k]]),exist_ok=True)

                write('{0}/D-{1}-{2}/POSCAR'.format(path,labels[oxygens[0][l]],labels[oxygens[1][k]]),H_lattice,sort=False)

    return traj, locations
