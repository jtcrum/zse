from zse.collections import framework
from zse.rings import get_osites
from zse.rings import get_tsites
import numpy as np

def label_osites(atoms, code):

    z = framework(code)
    osites,omult,first = get_osites(code)


    zcell = z.get_cell_lengths_and_angles()[:3]
    acell = atoms.get_cell_lengths_and_angles()[:3]
    repeat=[]
    for zc, ac in zip(zcell,acell):
        repeat.append(int(round(ac/zc)))

    z = z.repeat(repeat)
    oinds = [atom.index for atom in z if atom.symbol=='O']

    rp = np.prod(repeat)
    Dict = {}
    j=0
    for i in range(rp):
        for s,t in enumerate(osites):
            for q in range(omult[s]):
                Dict[oinds[j]]=t
                j+=1

    return Dict

def label_tsites(atoms, code):

    z = framework(code)
    tsites,tmult,first = get_tsites(code)
    tinds = [atom.index for atom in z if atom.symbol!='O']

    zcell = z.get_cell_lengths_and_angles()[:3]
    acell = atoms.get_cell_lengths_and_angles()[:3]
    repeat=[]
    for zc, ac in zip(zcell,acell):
        repeat.append(int(round(ac/zc)))
    z = z.repeat(repeat)
    tinds = [atom.index for atom in z if atom.symbol!='O']

    rp = np.prod(repeat)
    Dict = {}
    j=0
    for i in range(rp):
        for s,t in enumerate(tsites):
            for q in range(tmult[s]):
                Dict[tinds[j]]=t
                j+=1

    return Dict

def site_labels(atoms,code):
    '''
    This function will get the atom site labels (as defined by the IZA) for
    your atoms object. Be sure to remove any adsorbates from your zeolite
    framework before using this function or else it won't work. This function
    will work with T sites that have been exchanged for Al (or any atom).

    atoms: atoms object containing a zeolite you want labels for

    code: the zeolite framework code of your atoms object (i.e. 'CHA')
    '''


    tdict = label_tsites(atoms,code)
    odict = label_osites(atoms,code)
    all_labels = {**tdict,**odict}

    z = framework(code)
    zcell = z.get_cell_lengths_and_angles()[:3]
    acell = atoms.get_cell_lengths_and_angles()[:3]
    repeat=[]
    for zc, ac in zip(zcell,acell):
        repeat.append(int(round(ac/zc)))
    z = z.repeat(repeat)

    zo_inds = [atom.index for atom in z if atom.symbol =='O']
    zt_inds = [atom.index for atom in z if atom.symbol !='O']

    #z.set_tags(z_inds)
    poszo = z[zo_inds].get_scaled_positions()
    poszt = z[zt_inds].get_scaled_positions()

    Dict = {}
    for a in atoms:
        pa = a.scaled_position
        sym = a.symbol
        if sym == 'O':
            diffp = poszo - pa
            mags = []
            for d in diffp:
                mags.append(np.linalg.norm(d))
            ind = mags.index(min(mags))
            ind = zo_inds[ind]
            label = all_labels[ind]

        if sym != 'O':
            diffp = poszt - pa
            mags = []
            for d in diffp:
                mags.append(np.linalg.norm(d))
            ind = mags.index(min(mags))
            ind = zt_inds[ind]
            label = all_labels[ind]
        Dict[a.index]=label


    return Dict
