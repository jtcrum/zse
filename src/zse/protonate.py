from __future__ import annotations

from zse.proton_utilities import add_one_proton, add_two_protons, get_os_and_ts


def isolated(atoms, index, code, path=None):
    """
    Function to place a proton at each of the possible oxygens around a zeolite
    framework isolated T-site.

    INPUTS:
    atoms:      (ASE atoms object) Zeolite framework
    index:      (integer) index of the T-site to protonate
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')
    path:       (str) optional, if included structure files will be saved here

    OUTPUTS:
    traj:       (trajectory) Contains all the possible protonated structures.
    locations:  (list) Contains the location that correlates to each structure.
    """

    # first get all the oxygens and t sites needed
    oxygens, silicons = get_os_and_ts(atoms, index)

    # next add the proton to each oxygen
    traj, locations = add_one_proton(atoms, index, oxygens, silicons, code, path)

    return traj, locations


def paired(atoms, indices, code, path=None):
    """
    Function to enumerate all the possible proton locations at paired Al in a
    zeolite framework. This code won't work if the two Al are 1st nearest
    neighbor. Obey Löwenstein's rule for now I suppose.

    INPUTS:
    atoms:      (ASE atoms object) Zeolite framework
    indices:    (list of ints) indices of the T-sites to protonate
    code:       (str) IZA code for the zeolite you are using (i.e. 'CHA')
    path:       (str) optional, if included structure files will be saved here

    OUTPUTS:
    traj:       (trajectory) Contains all the possible protonated structures.
    locations:  (list) Contains the location that correlates to each structure.
    """

    # first get all the oxygens and silicons for each aluminum
    oxygens = []
    silicons = []
    for i in indices:
        tmpo, tmps = get_os_and_ts(atoms, i)
        oxygens.append(tmpo)
        silicons.append(tmps)

    # next add the protons to each possible oxygen, you will get 16 structures
    traj, locations = add_two_protons(atoms, indices, oxygens, silicons, code, path)

    return traj, locations
