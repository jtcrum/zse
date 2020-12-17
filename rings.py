'''
The goal of this module is to identify the rings associated with an oxygen or
tsite in any given zeolite framework. This method uses graph theory to find
neighbors of the specified atom, and then a depth first search for a cycle back
to that atom. To make the code as efficient as possible, its important to
include what types are rings are possible in that framework. This information is
stored within the collections module of this package. Check the examples page on
github (github.com/jtcrum/zse/examples) for specifics on how to use.
'''

import networkx as nx
from ase.io import read, write
from ase import neighborlist
import numpy as np
import math
from zse import substitute
from zse.collections import framework

def get_fwrings(code):
    from zse.collections import get_fwrings
    fw_rings = get_fwrings(code)
    return fw_rings

def get_orings_new(atoms,index,possible):
        cell = atoms.get_cell_lengths_and_angles()[:3]
        repeat = []
        possible = possible*2
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
        atoms2 = atoms.copy()
        atoms2 = atoms2.repeat(repeat)
        center = atoms2.get_center_of_mass()
        trans = center - atoms2.positions[index]
        atoms2.translate(trans)
        atoms2.wrap()
        atoms3 = atoms2.copy()
        cutoff = neighborlist.natural_cutoffs(atoms2,mult = 1.05)
        nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
        nl.update(atoms2)
        matrix = nl.get_connectivity_matrix(sparse=False)
        m = matrix.copy()
        G = nx.from_numpy_matrix(matrix)

        neighbs = nx.neighbors(G,index)
        fe = []
        for n in neighbs:
            if atoms2[n].symbol != 'O':
                fe.append(n)
        fe.append(index)
        rings = []
        for path in nx.all_simple_paths(G,index,fe[0],cutoff=max(possible)-1):
            rings.append(path)
        new_rings = []
        for r in rings:
            if len(r) in possible:
                new_rings.append(r)
        rings = new_rings

        delete = []
        for j,r in enumerate(rings):
            flag = False
            if len(r)>=12:
                for i in range(1,len(r)-3,2):
                    angle = atoms3.get_angle(r[i],r[i+2],r[i+4],mic=True)
                    if angle < 100:
                        delete.append(j)
                        break
        new_rings = []
        for j,r in enumerate(rings):
            if j not in delete:
                new_rings.append(r)
        rings = new_rings


        rings = remove_sec(rings)
        rings = remove_dups(rings)
        Class = []
        for r in rings:
            Class.append(int(len(r)/2))
        paths = rings
        paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
        Class.sort(reverse=True)

        keepers = []
        for i in paths:
            for j in i:
                if j not in keepers:
                    keepers.append(j)
        d = [atom.index for atom in atoms2 if atom.index not in keepers]
        del atoms2[d]


        return Class, paths, atoms2

def get_orings(atoms, index,possible):

    '''
    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    possible: (list) of the types of rings known to be present in the zeolite
              framework you are studying. This information is available on IZA
              or in the collections module of this package.
    Returns: Class - The size of rings associated with the oxygen.
             paths - The actual atom indices that compose those rings.
    '''

    cell = atoms.get_cell_lengths_and_angles()[:3]
    repeat = []
    possible = possible*2
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
    atoms2 = atoms.copy()
    atoms2 = atoms2.repeat(repeat)
    center = atoms2.get_center_of_mass()
    trans = center - atoms2.positions[index]
    atoms2.translate(trans)
    atoms2.wrap()
    atoms3 = atoms2.copy()
    cutoff = neighborlist.natural_cutoffs(atoms2,mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms2)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    fe = []
    for n in neighbs:
        if atoms2[n].symbol != 'O':
            fe.append(n)
    fe.append(index)

    tmpClass = []
    rings = []
    G2 = G.copy()
    G2.remove_edge(fe[0],fe[2])
    while len(tmpClass)<8:
        try:
            path = nx.shortest_path(G2,fe[0],fe[2])
        except:
            break
        length = len(path)
        if length in possible:
            tmpClass.append(int(len(path)/2))
            rings.append(path)
            if length == 18:
                G2.remove_edge(path[8],path[9])
            elif length < 16 and length > 6:
                if ((len(path)/2)%2)==0:
                    G2.remove_node(path[int(length/2-1)])
                elif ((len(path)/2)%2)!=0:
                    G2.remove_node(path[int(length/2)])
            elif length >=16:
                G2.remove_edge(path[int(len(path)/2-1)],path[int(len(path)/2)])
            if length == 8:
                G2.remove_node(path[4])
            if length == 6:
                G2.remove_node(path[3])
        else:
            if ((len(path)/2)%2)==0:
                G2.remove_node(path[int(length/2-1)])
            elif ((len(path)/2)%2)!=0:
                G2.remove_node(path[int(length/2)])

    tmpClass = []
    G2=G.copy()
    G2.remove_edge(fe[1],fe[2])
    rings2=[]
    while len(tmpClass)<8:
        try:
            path = nx.shortest_path(G2,fe[1],fe[2])
        except:
            break
        length = len(path)
        if length in possible:
            tmpClass.append(int(len(path)/2))
            rings2.append(path)
            if length == 18:
                G2.remove_edge(path[8],path[9])
            elif length < 16 and length > 6:
                if ((len(path)/2)%2)==0:
                    G2.remove_node(path[int(length/2-1)])
                elif ((len(path)/2)%2)!=0:
                    G2.remove_node(path[int(length/2)])
            elif length >=16:
                G2.remove_edge(path[int(len(path)/2-1)],path[int(len(path)/2)])
            if length == 8:
                G2.remove_node(path[4])
            if length == 6:
                G2.remove_node(path[3])
        else:
            if ((len(path)/2)%2)==0:
                G2.remove_node(path[int(length/2-1)])
            elif ((len(path)/2)%2)!=0:
                G2.remove_node(path[int(length/2)])

    rings = remove_sec(rings)
    rings2 = remove_sec(rings2)
    for i in rings2:
        rings.append(i)

    rings = remove_dups(rings)
    Class = []
    for r in rings:
        Class.append(int(len(r)/2))
    paths = rings
    paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
    Class.sort(reverse=True)

    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms2 if atom.index not in keepers]
    del atoms2[d]


    return Class, paths, atoms2

def get_rings(atoms, index):

    '''
    WARNING: This is old and does not work for all framework types.
    WARNING: Use the updated get_orings function below instead.


    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    '''

    cell = atoms.get_cell_lengths_and_angles()[:3]
    repeat = []

    for i,c in enumerate(cell):
        if c/2 <15:
            l = c
            re = 1
            while l/2 < 15:
                re +=1
                l = c*re

            repeat.append(re)
        else:
            repeat.append(1)
    atoms = atoms.repeat(repeat)
    center = atoms.get_center_of_mass()
    trans = center - atoms.positions[index]
    atoms.translate(trans)
    atoms.wrap()

    cutoff = neighborlist.natural_cutoffs(atoms,mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    for n in neighbs:
        if atoms[n].symbol == 'Si':
            fe = [n]
    fe.append(index)

    G.remove_edge(fe[0],fe[1])
    Class = []
    while len(Class)<6:
        try:
            path = nx.shortest_path(G,fe[0],fe[1])
        except:
            break
        Class.append(int(len(path)/2))
        for i in range(len(path)-3):
            G.remove_edge(path[i+1],path[i+2])
        Class.sort(reverse=True)
    return Class

def get_trings(atoms,index,possible):
    Class, paths, atoms2, repeat = tring_driver(atoms,index,possible)
    return Class, paths, atoms2

def get_tsites(code):
    from zse.collections import get_tsites
    z = framework(code)
    tsites,tmult = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol!='O']
    index = 0
    first_ts = []
    for i,m in enumerate(tmult):
        first_ts.append(tinds[index])
        index+=m
    return tsites,tmult,first_ts

def find_o_rings(G,index,possible):
    '''
    This is a helper function for the get_trings function.
    It won't do much on its own.
    '''
    rings = []
    neighbs = nx.neighbors(G,index)
    oxygen = []
    for n in neighbs:
        oxygen.append(n)
    for i in range(4):
        G2 = G.copy()
        G2.remove_edge(index, oxygen[i])
        tmp_class = []

        while len(tmp_class)<8:
            try:
                path = nx.shortest_path(G2,index,oxygen[i])
            except:
                break
            length = len(path)
            if length in possible:
                tmp_class.append(int(len(path)/2))
                rings.append(path)
                if length == 18:
                    G2.remove_edge(path[8],path[9])
                elif length < 16 and length > 6:
                    G2.remove_edge(path[3],path[4])
                elif length >=16:
                    G2.remove_edge(path[int(len(path)/2-1)],path[int(len(path)/2)])
                if length == 8:
                    G2.remove_node(path[3])
                if length == 6:
                    G2.remove_node(path[3])
            else:
                if (len(path)%2)==0:
                    G2.remove_node(path[int(length/2-1)])
                elif (len(path)%2)!=0:
                    G2.remove_node(path[int(length/2)])

    return rings

def remove_dups(rings):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    for i in range(len(rings)):
        for j in range((i+1), len(rings)):
            if i != j:
                st1 = set(rings[i])
                st2 = set(rings[j])
                if st1 == st2:
                    d.append(int(j))
    paths = []
    for i in range(len(rings)):
        if i not in d:
            paths.append(rings[i])
    return paths

def remove_sec(rings):
    '''
    This is a helper function for get_orings and get_trings.
    '''
    d = []
    count2 = np.zeros(len(rings))
    # for i in range(len(rings)):
    #     for j in range(i+1,len(rings)):
    #         if i!= j:
    #             ringi = rings[i]
    #             ringj = rings[j]
    #             ni = len(ringi)
    #             nj = len(ringj)
    #             if ni > nj and ni >= 16:
    #                 count=0
    #                 for rj in ringj:
    #                     if rj in ringi:
    #                         count+=1
    #                 if count > nj/2:
    #                     d.append(i)
    #             if nj >ni and nj >= 16:
    #                 count=0
    #                 for ri in ringi:
    #                     if ri in ringj:
    #                         count+=1
    #                 if count > ni/2:
    #                     d.append(j)
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            if i!= j:
                ringi = rings[i]
                ringj = rings[j]
                ni = len(ringi)
                nj = len(ringj)
                if ni > nj and ni >= 16 and nj > 6:
                    count=0
                    for rj in ringj:
                        if rj in ringi:
                            count+=1
                    if count == nj/2:
                        count2[i]+=1
                    elif count > nj/2:
                        count2[i]+=2
                if nj >ni and nj >= 16 and ni > 6:
                    count=0
                    for ri in ringi:
                        if ri in ringj:
                            count+=1
                    if count == ni/2:
                        count2[j]+=1
                    elif count > ni/2:
                        count2[j]+=2
                if ni > nj and nj in [6,8]:
                    count=0
                    for rj in ringj:
                        if rj in ringi:
                            count+=1
                    if count >= nj-2:
                        count2[i]+=2
                if nj > ni and ni in [6,8]:
                    count=0
                    for ri in ringi:
                        if ri in ringj:
                            count+=1
                    if count >= ni-2:
                        count2[j]+=2
    for i,c in enumerate(count2):
        if c >=2:
            d.append(i)

    paths = []
    for i in range(len(rings)):
        if i not in d:
            paths.append(rings[i])
    return paths

def tring_driver(atoms,index,possible,delete=True):
    '''
    atoms: ASE atoms object of the zeolite framework to be analyzed
    index: (integer) index of the atom that you want to classify
    possible: (list) of the types of rings known to be present in the zeolite
              framework you are studying. This information is available on IZA
              or in the collections module of this package.
    Returns: Class - The size of the rings associated with the desire T Site.
             Rings - The actual atom indices that compose those rings.
             atoms2 - An atoms object with the desired T Site changed to an
             Aluminum atom (just for visual purposes), and all atoms removed
             except for those that share a ring with the T Site provided.
    '''
    possible = possible*2
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
    rings = find_o_rings(G,index,possible)
    paths = remove_dups(rings)
    paths = remove_sec(paths)
    if delete == True:
        keepers = []
        for i in paths:
            for j in i:
                if j not in keepers:
                    keepers.append(j)
        d = [atom.index for atom in atoms2 if atom.index not in keepers]
        atoms2 = substitute.tsub(atoms2,index,'Al')
        del atoms2[d]

    Class = []
    for p in paths:
        Class.append(int(len(p)/2))

    paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
    Class.sort(reverse=True)

    return Class ,paths, atoms2, repeat

def unique_rings(code):
    z = framework(code)
    pr = get_fwrings(code)
    tsites,tmult, first = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol!='O']
    index = 0
    firstts = []
    for i,m in enumerate(tmult):
        firstts.append(tinds[index])
        index+=m
    allrings = []
    for f in firstts:
        c,r,ringatoms,repeat = tring_driver(z,f,pr,delete=False)
        for ring in r:
            allrings.append(ring)
    tinds = [atom.index for atom in ringatoms if atom.symbol!='O']
    rp = np.prod(repeat)
    Dict = {}
    j=0
    for i in range(rp):
        for s,t in enumerate(tsites):
            for q in range(tmult[s]):
                Dict[tinds[j]]=t
                j+=1
    ring_tsites = []
    for ring in allrings:
        tmp = []
        for i in ring:
            if ringatoms[i].symbol != 'O':
                tmp.append(Dict[i])
        ring_tsites.append(tmp)
    unique_tsites = {}
    unique_full = {}
    for i,r in enumerate(ring_tsites):
        length = len(r)
        if length not in unique_tsites:
            unique_tsites[length] = [r]
            unique_full[length] = [allrings[i]]
        else:
            unique_tsites[length].append(r)
            unique_full[length].append(allrings[i])
    trajectories = {}
    com = ringatoms.get_center_of_mass()
    for length in pr:
        ring_tlist = unique_tsites[length]
        ring_full  = unique_full[length]
        d = []
        for i in range(len(ring_tlist)):
            for j in range((i+1), len(ring_tlist)):
                st1 = ' '.join(map(str,ring_tlist[i]))
                st2 = ' '.join(map(str,ring_tlist[j]))
                st2_2 = ' '.join(map(str,reversed(ring_tlist[j])))
                if st2 in st1 + ' ' + st1 or st2_2 in st1 + ' ' + st1:
                    d.append(int(j))
        tmp1 = []
        tmp2 = []
        for i in range(len(ring_tlist)):
            if i not in d:
                tmp1.append(ring_tlist[i])
                tmp2.append(ring_full[i])
        unique_tsites[length] = tmp1
        unique_full[length] = tmp2
        traj = []
        for ring in tmp2:
            keepers = []
            atoms = ringatoms.copy()
            for i in ring:
                if i not in keepers:
                    keepers.append(i)
            d = [atom.index for atom in atoms if atom.index not in keepers]
            del atoms[d]
            position = atoms[0].position
            trans = com-position
            atoms.translate(trans)
            atoms.wrap()
            traj+=[atoms]
        trajectories[length]=traj


    return unique_tsites, unique_full, trajectories

'''
The following are developmental codes for testing
'''

def test_orings(atoms,index,possible):
    cell = atoms.get_cell_lengths_and_angles()[:3]
    repeat = []
    possible = possible*2
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
    atoms2 = atoms.copy()
    atoms2 = atoms2.repeat(repeat)
    center = atoms2.get_center_of_mass()
    trans = center - atoms2.positions[index]
    atoms2.translate(trans)
    atoms2.wrap()

    cutoff = neighborlist.natural_cutoffs(atoms2,mult = 1.05)
    nl = neighborlist.NeighborList(cutoffs = cutoff, self_interaction=False, bothways = True)
    nl.update(atoms2)
    matrix = nl.get_connectivity_matrix(sparse=False)
    m = matrix.copy()
    G = nx.from_numpy_matrix(matrix)

    neighbs = nx.neighbors(G,index)
    fe = []
    for n in neighbs:
        if atoms2[n].symbol != 'O':
            fe.append(n)
    fe.append(index)

    tmpClass = []
    rings = []
    G2 = G.copy()
    G2.remove_edge(fe[0],fe[2])
    pr = sorted(possible)
    rings = []

    for q in pr:
        tmprings = []
        paths = nx.all_simple_paths(G2,fe[0],fe[2],cutoff=q-1)
        for p in paths:
            tmprings.append(p)
        for r in tmprings:
            length =len(r)
            if length in possible:
                if ((len(r)/2)%2)==0:
                    try:
                        G2.remove_node(r[int(length/2-1)])
                    except:
                        2+2
                elif ((len(r)/2)%2)!=0:
                    try:
                        G2.remove_node(r[int(length/2)])
                    except:
                        2+2

                rings.append(r)
    rings = remove_dups(rings)
    rings = remove_sec(rings)
    Class = []
    for r in rings:
        Class.append(int(len(r)/2))
    paths = rings
    paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
    Class.sort(reverse=True)

    keepers = []
    for i in paths:
        for j in i:
            if j not in keepers:
                keepers.append(j)
    d = [atom.index for atom in atoms2 if atom.index not in keepers]
    del atoms2[d]


    return Class, paths, atoms2

def get_all_rings(code):
    '''
    For developmental testing only.
    '''
    z = framework(code)
    pr = get_fwrings(code)
    tsites,tmult, first = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol!='O']
    index = 0
    firstts = []
    for i,m in enumerate(tmult):
        firstts.append(tinds[index])
        index+=m
    allrings = []
    for f in firstts:
        c,r,ringatoms,repeat = all_trings(z,f,pr,delete=False)
        for ring in r:
            allrings.append(ring)
    tinds = [atom.index for atom in ringatoms if atom.symbol!='O']
    rp = np.prod(repeat)
    Dict = {}
    j=0
    for i in range(rp):
        for s,t in enumerate(tsites):
            for q in range(tmult[s]):
                Dict[tinds[j]]=t
                j+=1
    ring_tsites = []
    for ring in allrings:
        tmp = []
        for i in ring:
            if ringatoms[i].symbol != 'O':
                tmp.append(Dict[i])
        ring_tsites.append(tmp)
    unique_tsites = {}
    unique_full = {}
    for i,r in enumerate(ring_tsites):
        length = len(r)
        if length not in unique_tsites:
            unique_tsites[length] = [r]
            unique_full[length] = [allrings[i]]
        else:
            unique_tsites[length].append(r)
            unique_full[length].append(allrings[i])
    trajectories = {}
    com = ringatoms.get_center_of_mass()
    for length in pr:
        ring_tlist = unique_tsites[length]
        ring_full  = unique_full[length]
        d = []
        for i in range(len(ring_tlist)):
            for j in range((i+1), len(ring_tlist)):
                st1 = ' '.join(map(str,ring_tlist[i]))
                st2 = ' '.join(map(str,ring_tlist[j]))
                st2_2 = ' '.join(map(str,reversed(ring_tlist[j])))
                if st2 in st1 + ' ' + st1 or st2_2 in st1 + ' ' + st1:
                    d.append(int(j))
        tmp1 = []
        tmp2 = []
        for i in range(len(ring_tlist)):
            if i not in d:
                tmp1.append(ring_tlist[i])
                tmp2.append(ring_full[i])
        unique_tsites[length] = tmp1
        unique_full[length] = tmp2
        traj = []
        for ring in tmp2:
            keepers = []
            atoms = ringatoms.copy()
            for i in ring:
                if i not in keepers:
                    keepers.append(i)
            d = [atom.index for atom in atoms if atom.index not in keepers]
            del atoms[d]
            position = atoms[0].position
            trans = com-position
            atoms.translate(trans)
            atoms.wrap()
            traj+=[atoms]
        trajectories[length]=traj


    return unique_tsites, unique_full, trajectories

def all_trings(atoms,index,possible,delete=True):
    '''
    For developmental testing purposes only.
    '''
    possible = possible*2
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
    rings = find_simple_o_rings(G,index,possible)
    paths = remove_dups(rings)
    # paths = remove_sec(paths)
    if delete == True:
        keepers = []
        for i in paths:
            for j in i:
                if j not in keepers:
                    keepers.append(j)
        d = [atom.index for atom in atoms2 if atom.index not in keepers]
        atoms2 = substitute.tsub(atoms2,index,'Al')
        del atoms2[d]

    Class = []
    for p in paths:
        Class.append(int(len(p)/2))

    paths = [x for _,x in sorted(zip(Class,paths),reverse=True)]
    Class.sort(reverse=True)

    return Class ,paths, atoms2, repeat

def find_simple_o_rings(G,index,possible):
    '''
    This is a helper function for the get_trings function.
    It won't do much on its own.
    '''
    rings = []
    neighbs = nx.neighbors(G,index)
    oxygen = []
    for n in neighbs:
        oxygen.append(n)
    for i in range(4):
        G2 = G.copy()
        G2.remove_edge(index, oxygen[i])
        tmp_class = []

        for path in nx.all_simple_paths(G2,index,oxygen[i],cutoff=max(possible)-1):
            rings.append(path)

    return rings

def get_unique_trings(code,ring_size):
    z = framework(code)
    pr = get_fwrings(code)
    tsites,tmult = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol!='O']
    index = 0
    firstts = []
    for i,m in enumerate(tmult):
        firstts.append(tinds[index])
        index+=m
    allrings = []
    for f in firstts:
        c,r,ringatoms,repeat = tring_driver(z,f,pr,delete=False)
        for ring in r:
            allrings.append(ring)
    tinds = [atom.index for atom in ringatoms if atom.symbol!='O']
    rp = np.prod(repeat)
    Dict = {}
    j=0
    for i in range(rp):
        for s,t in enumerate(tsites):
            for q in range(tmult[s]):
                Dict[tinds[j]]=t
                j+=1
    ring_tsites = []
    for ring in allrings:
        tmp = []
        for i in ring:
            if ringatoms[i].symbol != 'O':
                tmp.append(Dict[i])
        ring_tsites.append(tmp)
    desired_rings_tsites = []
    desired_rings_full = []
    for i,r in enumerate(ring_tsites):
        if len(r)==ring_size:
            desired_rings_tsites.append(r)
            desired_rings_full.append(allrings[i])
    unique_tsites = []
    unique_full = []
    d = []
    for i in range(len(desired_rings_tsites)):
        for j in range((i+1), len(desired_rings_tsites)):
            if i != j:
                st1 = ' '.join(map(str,desired_rings_tsites[i]))
                st2 = ' '.join(map(str,desired_rings_tsites[j]))
                st2_2 = ' '.join(map(str,reversed(desired_rings_tsites[j])))
                # st1 = set(desired_rings_tsites[i])
                # st2 = set(desired_rings_tsites[j])
                # if st1 == st2:
                if st2 in st1 + ' ' + st1 or st2_2 in st1 + ' ' + st1:
                    d.append(int(j))
    for i in range(len(desired_rings_tsites)):
        if i not in d:
            unique_tsites.append(desired_rings_tsites[i])
            unique_full.append(desired_rings_full[i])

    traj = []
    com = ringatoms.get_center_of_mass()
    for ring in unique_full:
        keepers = []
        atoms = ringatoms.copy()
        for i in ring:
            if i not in keepers:
                keepers.append(i)
        d = [atom.index for atom in atoms if atom.index not in keepers]
        del atoms[d]
        position = atoms[0].position
        trans = com-position
        atoms.translate(trans)
        atoms.wrap()
        traj+=[atoms]

    return traj, unique_tsites
