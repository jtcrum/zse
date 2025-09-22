from zse.cation import divalent, monovalent
from zse.collections import framework
from zse.substitute import tsub

CHA_SINGLE_IDX = 101
CHA_DOUBLE_IDX = [98, 101]
MONOVALENT_ION = "Na"
DIVALENT_ION = "Cu"


def test_divalent():
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute two Si atoms with Al
    atoms = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Place a divalent cation
    traj = divalent(atoms, DIVALENT_ION)

    # Check that 12 structures are generated
    assert len(traj) == 12


def test_monovalent_all_rings():
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Place a monovalent cation
    traj, locations = monovalent(atoms, CHA_SINGLE_IDX, MONOVALENT_ION)

    assert len(traj) == 4
    assert len(locations) == 4


def test_monovalent_specific_rings():
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Place a monovalent cation
    traj, locations = monovalent(atoms, CHA_SINGLE_IDX, MONOVALENT_ION, included_rings=[8, 6, 4])

    assert len(traj) == 6
    assert len(locations) == 6
