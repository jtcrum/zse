__all__ = ['read_cif','cif_site_labels']

from ase.io import read
from ase.spacegroup import spacegroup
import sys
import os
import logging
from math import *
import numpy as np
import pkg_resources
import warnings
warnings.filterwarnings("ignore")


path = '.temp_files/'
filepath = pkg_resources.resource_filename(__name__,path)

'''
NOTE ABOUT CIF FILE FORMATS:
CIFs must include '_symmetry_Int_Taables_number' to be read by ASE.
If this is not included please edit your CIF file to include this information.
'''

def get_atom_lines(alllines):
    order = []
    for i,line in enumerate(alllines):
        if '_atom' in line:
            order.append(line)
            start = i+1
    end = None
    for i,line in enumerate(alllines[start:]):

        if len(line.split()) == 0:
            end = start+i-1
            break
    if not end:
        end = len(alllines)-1

    new_order = []
    for i,o in enumerate(order):
        if 'site_label' in o:
            new_order.append(i)
        if 'site_type_symbol' in o:
            new_order.append(i)
        if 'fract_x' in o:
            new_order.append(i)
        if 'fract_y' in o:
            new_order.append(i)
        if 'fract_z' in o:
            new_order.append(i)

    return start,end,new_order

def fix_cif(cif):
    f = open(cif,"r")
    alllines = f.readlines()
    f.close()

    for i, line in enumerate(alllines):
        if 'IT_coordinate_system_code' in line:
            fields = line.split()
            alllines[i] = '_symmetry_space_group_setting {0} \n'.format(fields[-1])

        if '_atom_site_type_symbol' in line and '_atom_site_label' in alllines[i+1]:
            alllines[i],alllines[i+1] = alllines[i+1],alllines[i]

    file_name = cif.rstrip('.cif')
    temp_file = '{0}/{1}_temp.cif'.format(filepath,file_name.split('/')[-1])
    f = open(temp_file,"w")
    f.writelines(alllines)
    f.close()
    atoms = read(temp_file);
    os.remove(temp_file)
    return atoms, alllines

def get_tsites(cif):
    from ase.geometry import get_distances
    tsites = []
    tpos = []
    z,alllines = fix_cif(cif)
    si = [atom.index for atom in z if atom.symbol!='O']
    start,end,order = get_atom_lines(alllines)
    for line in alllines[start:end+1]:
        if 'Si' in line or 'T' in line:
            line = line.split()
            temp_label = line[order[0]]
            if not any(str.isdigit(c) for c in temp_label):
                temp_label = line[order[1]]
            if 'Si' in temp_label:
                temp_label = temp_label.replace('Si','T')
            tsites.append(temp_label)
            pos = [float(line[order[2]]),float(line[order[3]]),float(line[order[4]])]
            tpos.append([round(num,2) for num in pos])

    tpos = np.array(tpos)
    pos = z[si].get_scaled_positions()
    tinds = []
    tmults = []
    t_class = []
    for tp in tpos:
        for i,p in enumerate(pos):
            p = [round(num,2) for num in p]
            diff = abs(tp-p)
            if sum(diff) <= 0.03:
                tinds.append(si[i])

    for i in range(1,len(tsites)):
        tmults.append(tinds[i]-tinds[i-1])
    tmults.append(si[-1]-tinds[-1]+1)

    #
    # si = [atom.index for atom in z if atom.symbol=='Si']
    # o  = [atom.index for atom in z if atom.symbol=='O']
    # si_pos = z[si].positions
    # cell = z.cell
    # distances = get_distances(si_pos,si_pos,cell=cell,pbc=[1,1,1])[1]
    #
    # for i in tinds:
    #     orig_ind = si.index(i)
    #     dists = sorted(distances[orig_ind])
    #     t_class.append([round(num,2) for num in dists])
    #
    #
    # for i,d in enumerate(t_class):
    #     for j,t in enumerate(distances):
    #         dist = [round(num,2) for num in sorted(t)]
    #         if np.array_equal(dist,d):
    #             dist = [round(num,2) for num in sorted(t)]
    #             d = np.array(d)
    #             dist = np.array(dist)
    #             diff = abs(d - dist)
    #             if sum(diff) <= 0.1:
    #                 tmults[i]+=1

    n = len(si)
    sn = sum(tmults)
    if n != sn:
        print('Something Went Wrong With T Sites')
    return tsites, tmults, tinds

def get_osites(cif):
    from ase.geometry import get_distances
    osites = []
    opos = []
    z,alllines = fix_cif(cif)
    start,end,order = get_atom_lines(alllines)
    for line in alllines[start:end+1]:
        if 'O' in line:
            line = line.split()
            temp_label = line[order[0]]
            if not any(str.isdigit(c) for c in temp_label):
                temp_label = line[order[1]]
            osites.append(temp_label)
            pos = [float(line[order[2]]),float(line[order[3]]),float(line[order[4]])]
            opos.append([round(num,2) for num in pos])
    opos = np.array(opos)
    pos = z.get_scaled_positions()
    oinds = []
    omults = []
    o_class = []

    si = [atom.index for atom in z if atom.symbol=='Si']
    o  = [atom.index for atom in z if atom.symbol=='O']
    o_pos = z[o].get_scaled_positions()
    for op in opos:
        for i,p in enumerate(o_pos):
            p = np.array([round(num,2) for num in p])
            diff = abs(op-p)
            if sum(diff) <= 0.02:
                oinds.append(o[i])

    for i in range(1,len(osites)):
        omults.append(oinds[i]-oinds[i-1])
    omults.append(o[-1]-oinds[-1]+1)

    # all_pos = z.positions
    # o_pos = z[o].positions
    # si_pos = z[si].positions
    # cell = z.cell
    # distances = get_distances(o_pos,all_pos,cell=cell,pbc=[1,1,1])[1]
    #
    # for i in oinds:
    #     orig_ind = o.index(i)
    #     dists = sorted(distances[orig_ind])
    #     o_class.append([round(num,2) for num in dists])
    #
    # for i,d in enumerate(o_class):
    #     for j,t in enumerate(distances):
    #         dist = [round(num,2) for num in sorted(t)]
    #         d = np.array(d)
    #         dist = np.array(dist)
    #         diff = abs(d - dist)
    #         if sum(diff) <= 0.05:
    #             omults[i]+=1

    n = len(o)
    sn = sum(omults)
    if n != sn:
        print('Something Went Wrong With O Sites')
    return osites, omults, oinds

def read_cif(cif):
    atoms, alllines = fix_cif(cif)
    ts,tm,tinds = get_tsites(cif)
    os,om,oinds = get_osites(cif)
    return atoms,ts,tm,tinds,os,om,oinds

def cif_site_labels(cif):
    atoms,ts,tm,tinds,os,om,oinds = read_cif(cif)
    labels = {}
    for i,t in enumerate(ts):
        for j in range(tm[i]):
            labels[tinds[i]+j] = t

    for i,o in enumerate(os):
        for j in range(om[i]):
            labels[oinds[i]+j] = o

    return labels

''' DEPRECRATED FUNCTIONS'''

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
    f = open(cif,"r")
    alllines = f.read()
    f.close()

    for i, line in enumerate(alllines):
        if 'IT_coordinate_system_code' in line:
            fields = line.split()
            alllines[i] = '_symmetry_space_group_setting {0}'.format(fields[-1])


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
