from ase import Atoms

from zse.proton_utilities import add_one_proton, add_two_protons, get_os_and_ts

__all__ = ["isolated", "paired"]


def isolated(
    atoms: Atoms, index: int, code: str, path: str | None = None
) -> tuple[list[Atoms], list[str]]:
    """Place a proton at each of the possible oxygens around a zeolite
    framework isolated T-site.

    Args:
        atoms (Atoms): The ASE Atoms object representing the structure.
        index (int): The index of the T-site to protonate.
        code (str): The IZA code for the zeolite you are using (i.e. 'CHA').
        path (str | None, optional): The directory path to save the modified structures.
            Defaults to None.

    Returns:
        tuple[list[Atoms], list[str]]: A list of modified Atoms objects and a list of
            location labels.
    """

    # first get all the oxygens and t sites needed
    oxygens, silicons = get_os_and_ts(atoms, index)

    # next add the proton to each oxygen
    traj, locations = add_one_proton(atoms, index, oxygens, silicons, code, path)

    return traj, locations


def paired(
    atoms: Atoms, indices: list[int], code: str, path: str | None = None
) -> tuple[list[Atoms], list[str]]:
    """Enumerate all the possible proton locations at paired Al in a
    zeolite framework. This code won't work if the two Al are 1st nearest
    neighbor. Obey Löwenstein's rule for now I suppose.

    Args:
        atoms (Atoms): The ASE Atoms object of the zeolite structure.
        indices (list[int]): The indices of the T-sites to protonate.
        code (str): The IZA code for the zeolite you are using (e.g., 'CHA').
        path (str | None, optional): The directory path to save the modified structures.
            Defaults to None.

    Returns:
        tuple[list[Atoms], list[str]]: A list of modified Atoms objects and a list of
            location labels.
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
