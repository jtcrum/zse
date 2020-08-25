from ase.io import read, write
from ase import Atoms, Atoms
import math
import os.path
import numpy as np
import os
from ase.build import molecule

def monovalent(atoms,symbol,path=None):
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

return traj

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
