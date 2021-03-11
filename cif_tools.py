__all__ = ['get_indices']

from ase.io import read
import sys
import logging
from math import *
import numpy as np

def float_with_error(x):
    """
    some value in cif accompanies error like "1.234(5)
    """
    if "?" in x:
        return 0
    pos = x.find("(")
    if pos >= 0:
        x = x[:pos]
    return float(x)

def get_mults(cif):

    # read the cif file

    F = open(cif,"r")
    alllines = F.readlines()
    F.close()

    # Parse out data from the cif file

    for i,line in enumerate(alllines):
        if '_cell_length_a' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            La = field
        if '_cell_length_b' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            Lb = field
        if '_cell_length_c' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            Lc = field
        if '_cell_angle_alpha' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            alpha = field
        if '_cell_angle_beta' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            beta = field
        if '_cell_angle_gamma' in line:
            fields = line.split()
            field = fields[-1]
            field = float_with_error(field)
            gamma = field
        if '_space_group_symop' in line or '_symmetry_equiv_pos' in line or '_space_group' in line:
            n = i

    lastline = len(alllines)

    loops = []
    for i,line in enumerate(alllines):
        if 'loop' in line:
            loops.append(i)

    ops = []
    for i in range(n+1,loops[1]):
        n+=1
        line = alllines[i]
        if 'x' in line or 'X' in line:
            ops.append(line.replace("'",''))

    for i in range(len(ops)):
        ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
        ops[i] = ops[i].replace("1/", "1./")
        ops[i] = ops[i].replace("2/", "2./")
        ops[i] = ops[i].replace("3/", "3./")
        ops[i] = ops[i].replace("4/", "4./")
        ops[i] = ops[i].replace("5/", "5./")
        ops[i] = ops[i].replace("6/", "6./")
        ops[i] = ops[i].replace("7/", "7./")
        ops[i] = ops[i].replace("8/", "8./")
        ops[i] = ops[i].replace("9/", "9./")

    osites = []
    tsites = []
    atoms = []
    for j in range(n,lastline):
        line = alllines[j]
        if '_' not in line:
            fields = line.split()
            if len(fields) >3:
                tmp = (fields[0],float(fields[2]),float(fields[3]),float(fields[4]))
                if 'O' in fields[0]:
                    osites.append(fields[0])
                if 'T' in fields[0]:
                    tsites.append(fields[0])
                atoms.append(tmp)
    for i in range(len(atoms)):
        (name,xn,yn,zn) = atoms[i]
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0
        atoms[i] = (name,xn,yn,zn)

    # perfrom symmetry operations

    label_list = []
    symbols = []
    positions = []

    for i in atoms:
        label_list.append(i[0])
    eps = 0.01
    imax = len(atoms)
    i=0
    while (i<imax):
        label,x,y,z=atoms[i]
        for op in ops:
            op = op.replace("'",'')
            op = op.lower()
            xn,yn,zn = eval(op)

            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0

            new_atom = True
            for at in atoms:
                if (abs(at[1]-xn) < eps and abs(at[2]-yn) < eps and abs(at[3]-zn) < eps):
                    new_atom = False
                if new_atom:
                    p1 = np.array([at[1],at[2],at[3]])
                    p2 = np.array([xn,yn,zn])
                    diff = abs(p1-p2)
                    diff = np.round(diff,2)
                    count = np.count_nonzero(diff)
                    if count ==1 and 1 in diff:
                        new_atom = False
            if new_atom:

                atoms.append( (label,xn,yn,zn) )
                label_list.append(label)
        i += 1
        imax =len(atoms)
    #atoms2 = Atoms(symbols,scaled_positions=positions,cell = [La,Lb,Lc,alpha,beta,gamma])
    # count up the osits
    label_list = sorted(label_list)
    omults = []
    for o in osites:
        count = label_list.count(o)
        omults.append(count)
    tmults = []
    for t in tsites:
        count = label_list.count(t)
        tmults.append(count)
    return tsites, tmults, osites, omults

def get_indices(cif):
    '''
    This is a tool that will read a CIF file and return the unique T-sites,
    their multiplicities, and an example atom index.

    It also does the same for the unique O-sites in the framework.

    This tool only works on CIFs that are formatted the same way as the IZA
    Structure Database CIFs.
    '''

    tsites, tmults, osites, omults = get_mults(cif)

    atoms = read(cif)

    oinds = [atom.index for atom in atoms if atom.symbol=='O']
    index  = 0
    first_os = []
    for i,m in enumerate(omults):
        first_os.append(oinds[index])
        index+=m

    tinds = [atom.index for atom in atoms if atom.symbol !='O']
    index = 0
    first_ts = []
    for i,m, in enumerate(tmults):
        first_ts.append(tinds[index])
        index+=m
    return tsites,tmults,first_ts, osites, omults, first_os
