__all__ = ['divalent','monovalent']

from zse.cation_utilities import *
from zse.ring_utilities import *
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

def monovalent(atoms,index,symbol,code,included_rings=None,path=None,bvect=None):

    '''
    This code has been updated to place the ion inside each of the rings
    associated with the T site. The rings are found using the rings module of
    ZSE.

    INPUTS:
    atoms = ASE atoms object of the zeolite framework
    index = Index of the t site that the cation will be associated with (int)
    symbol = Elemental symbol of the cation you want to use, i.e. 'Na' (str)
    code = Framework code of the zeolite you are using, i.e. 'CHA' (str)
    included_rings (optional) = List of ints for rings you want to include
                                if not specified all rings larger than 4-MR will
                                be inlcuded.
    path (optional) = Path for which you would like the structure files saved.
                      If not included, structure files will not be saved.
    bvect (optional) = Manually specify the bond length between the cation and 
                        atom index

    OUTPUTS:
    traj = ASE trajectory of all the structures generated. You can view traj
           with ase.visualize.view.
    locations = List of all the rings that the ion was placed in. Correlates to
                the images in the trajectory.
    '''

    # let's import the modules we will need to make this work
    from zse.collections import get_ring_sizes
    import networkx as nx
    from ase.data import covalent_radii, chemical_symbols

    # I will use these radii to approximate bond lengths
    radii = {chemical_symbols[i]: covalent_radii[i] for i in range(len(chemical_symbols))}

    # get all the rings associated with the T site
    # this follows the same steps as rings.get_trings()
    ring_sizes = get_ring_sizes(code)*2
    max_ring = max(ring_sizes)
    G, large_atoms, repeat = atoms_to_graph(atoms,index,max_ring)
    import networkx as nx
    paths = []
    for n in nx.neighbors(G,index):
        paths = paths+get_paths(G,n,ring_sizes)
    paths = remove_non_rings(large_atoms, paths)

    # which rings should be included
    if included_rings == None:
        included_rings = []
        for p in ring_sizes:
            if p > 8:
                included_rings.append(p)
    else:
        included_rings = [x*2 for x in included_rings]

    # get a list of all the rings and their sizes present
    Class, class_count, paths = count_rings(paths)

    # add the cation to each ring, put structure in a trajectory
    traj, locations = add_cation(atoms,large_atoms,radii,index,symbol,paths,included_rings,class_count,path,bvect)

    return traj, locations
