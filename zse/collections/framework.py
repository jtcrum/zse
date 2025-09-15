import importlib.resources as pkg_resources

from ase import Atoms
from ase.db import connect

path = "frameworks.db"
filepath = pkg_resources.files(__name__).joinpath(path)


def framework_db(code: str) -> Atoms:
    """Get the atoms object for a given zeolite framework code from the database.

    Args:
        code (str): Zeolite framework code.

    Returns:
        Atoms: Atoms object containing the zeolite framework.
    """
    db = connect(filepath)
    atoms = db.get_atoms(fw=code)
    return atoms


def get_ring_sizes_db(code: str) -> list[int]:
    """Get the ring sizes for a given zeolite framework code from the database.

    Args:
        code (str): Zeolite framework code.

    Returns:
        list[int]: A list of ring sizes for the specified framework code.
    """
    db = connect(filepath)
    fw_rings = db.get(fw=code).data.rings
    return fw_rings


def get_tsites_db(code: str) -> tuple[list[str], list[int], list[int]]:
    """Get the T sites for a given zeolite framework code from the database.

    Args:
        code (str): Zeolite framework code.

    Returns:
        tuple[list[str], list[int], list[int]]: A tuple containing a list of T site labels,
            a list of their multiplicities, and a list of the first T atom indices for each site.
    """
    db = connect(filepath)
    tsites = db.get(fw=code).data.tsites
    tmult = db.get(fw=code).data.tmult
    return tsites, tmult


def get_osites_db(code: str) -> tuple[list[str], list[int], list[int]]:
    """
    Get the oxygen sites for a given zeolite framework code from the database.
    Args:
        code (str): Zeolite framework code.

    Returns:
        tuple[list[str], list[int]]: A tuple containing a list of oxygen site labels and a
            list of their multiplicities.
    """
    db = connect(filepath)
    osites = db.get(fw=code).data.osites
    omult = db.get(fw=code).data.omult
    return osites, omult


def get_all_fws_db() -> list[str]:
    """Get a list of all zeolite framework codes in the database.

    Returns:
        list[str]: A list of all zeolite framework codes in the database.
    """
    db = connect(filepath)
    fws = []
    for row in db.select():
        fws.append(row.fw)  # noqa: PERF401
    return fws
