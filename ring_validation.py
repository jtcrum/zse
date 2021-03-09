__all__ = ['sp','d2','sphere','cross_distance']
from zse.ring_utilities import *

def sp(G,paths):
    '''
    Custum method for determing valid rings by ensuring for paths of 8 T-sites
    or larger, the current path is the shortest possible way to connect each
    node long the path. Alternate paths must be shorter than the current path
    to be considered invalid.

    7 T-site paths and smaller are handled slightly
    different where if the alrternate path is equal in length to the current
    path it is considered non valid.
    '''

    import networkx as nx
    valid_paths = []
    for p in paths:
        l = len(p)
        flag = True
        if l/2 < 8:
            for j in range(1,l-3,2):
                for k in range(j+4,l,2):
                    shortest_path = nx.shortest_path(G,p[j],p[k])
                    for r in shortest_path:
                        if r not in p:
                            lengths = [k-j+1,l-k+j+1]
                            if len(shortest_path)<min(lengths):
                                flag = False
                                break
        elif l/2 >=8:
            for j in range(1,l-3,2):
                for k in range(j+2,l,2):
                    shortest_path = nx.shortest_path(G,p[j],p[k])
                    for r in shortest_path:
                        if r not in p:
                            lengths = [k-j+1,l-k+j+1]
                            if len(shortest_path)<max(lengths):
                                flag = False
                                break
            # if not flag:
            #     break
        if flag:
            valid_paths.append(p)

    return valid_paths

def sphere(atoms,paths,cutoff):
    '''
    Method for determining valid rings by ensuring that non ring atoms are not
    within a cutoff distance of the center of mass of the path.
    '''
    from ase.geometry import get_distances
    valid_paths = []
    for p in paths:
        flag = True
        if len(p) > 10:
            ring_atoms = atoms[p]
            com = ring_atoms.get_center_of_mass()
            positions = atoms.get_positions()
            distances = get_distances(com,positions)[1][0]

            for i,d in enumerate(distances):
                if d < cutoff and i not in p:
                    flag = False
                    break
        if flag:
            valid_paths.append(p)

    return valid_paths

def cross_distance(atoms, paths):
    '''
    This validation method uses cross Si-Si distances of the path to determine
    if the path is a valid ring or not. The cutoffs for the cross Si-Si
    distances were determined by analyzing many valid and invalid rings.
    This method is biased by human interpretation of what constitutes a ring.
    '''
    import math

    delete = []
    for j,r in enumerate(paths):
        n = int(len(r)/2)
        cutoff = n/4
        if cutoff < 2:
            cutoff = 2
        if n%2 == 0 and n > 5:
            distances = []
            inner_flag = False
            for x in range(1,n,2):
                dist = atoms.get_distance(r[x],r[x+n],mic=True)
                distances.append(dist)
                if dist < n-cutoff:
                    delete.append(j)
                    inner_flag = True
                    break
            if inner_flag == False:
                outer_flag = False
                for d in distances:
                    if d > n - math.floor(n/6):
                        outer_flag = True
                        break
                if outer_flag == False:
                    delete.append(j)
        if n%2 != 0 and n > 5:
            r2 = r.copy()
            r2.append(r[:2])
            for x in range(1,n,2):
                dist1 = atoms.get_distance(r2[x],r2[x+n-1],mic=True)
                dist2 = atoms.get_distance(r2[x],r2[x+n+1],mic=True)
                if dist1 < n-cutoff or dist2 < n-cutoff:
                    delete.append(j)
                    break

    tmp_paths = paths.copy()
    paths = []
    for j,r in enumerate(tmp_paths):
        if j not in delete:
            paths.append(r)

    return paths

def d2(G,paths):
    '''
    Method for determinging valid rings presented by Goetzke, K.; Klein, H.-J.
    (DOI: 10.1016/0022-3093(91)90145-V) and implemented by
    Sastre, G; Corma, A. (DOI: 10.1021/jp8100128)
    A valid ring is a path that cannote be decomposed into two smaller rings.
    Results found with this method match results from the Sastre & Corma paper.
    '''
    import networkx as nx
    valid_paths = []
    for path in paths:
        FLAG = False
        path2 = path + path
        l = len(path)
        if l > 8:
            for j in range(1,l,2):
                for k in range(j+4,int(l/2)+j+1,2):
                    p1 = path2[j:k+1]
                    p2 = path2[k:l+j+1]
                    l1 = len(p1)
                    l2 = len(p2)
                    if l1 < l2:
                        G2 = G.copy()
                        for x in p1[1:-1]:
                            G2.remove_node(x)
                        sp = nx.shortest_path(G2,p1[0],p1[-1])
                        con1 = p2 + sp[1:-1]
                        if len(con1) < l:
                            FLAG = True
                        if len(con1) == l:
                            FLAG,q = is_valid(G,con1)
                    if l2 < l2:
                        G2 = G.copy()
                        for x in p2[1:-1]:
                            G2.remove_node(x)
                        sp = nx.shortest_path(G2,p2[0],p2[-1])
                        con1 = p1 + sp[1:-1]
                        if len(con1) < l:
                            FLAG = True
                        if len(con1) == l:
                            FLAG,q = is_valid(G,con1)
                    if FLAG:
                        break
                if FLAG:
                    break
        if not FLAG:
            valid_paths.append(path)
    return valid_paths
