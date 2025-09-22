import numpy as np

from zse.collections import framework
from zse.proton_utilities import add_one_proton, add_two_protons, get_os_and_ts
from zse.substitute import tsub

CHA_SINGLE_IDX = 101
CHA_DOUBLE_IDX = [98, 101]


def test_get_os_and_ts():
    """Test getting oxygen and silicon atoms around a T-site."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Get oxygens and silicons around a T-site
    oxygens, silicons = get_os_and_ts(atoms, CHA_SINGLE_IDX)

    # Check that arrays are returned
    assert isinstance(oxygens, np.ndarray)
    assert isinstance(silicons, np.ndarray)

    # Check that we get the expected number of oxygens (4 for a T-site)
    assert len(oxygens) == 4

    # Check that we get some silicons
    assert len(silicons) > 0

    # Check that all returned indices are valid
    assert all(0 <= idx < len(atoms) for idx in oxygens)
    assert all(0 <= idx < len(atoms) for idx in silicons)

    # Check that oxygens are actually oxygen atoms
    for oxy_idx in oxygens:
        assert atoms[oxy_idx].symbol == "O"

    # Check that silicons are actually silicon atoms
    for si_idx in silicons:
        assert atoms[si_idx].symbol == "Si"


def test_add_one_proton():
    """Test adding a single proton to the structure."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Get oxygens and silicons around the substituted site
    oxygens, silicons = get_os_and_ts(atoms, CHA_SINGLE_IDX)

    # Add one proton
    traj, locations = add_one_proton(atoms, CHA_SINGLE_IDX, oxygens, silicons, "CHA")

    # Check that structures and locations are returned
    assert isinstance(traj, list)
    assert isinstance(locations, list)
    assert len(traj) == len(locations)

    # Should get 4 structures (one for each oxygen)
    assert len(traj) == 4

    # Check that each structure has a proton added
    for atoms_with_proton in traj:
        assert "H" in atoms_with_proton.get_chemical_symbols()
        assert len(atoms_with_proton) == len(atoms) + 1
        # Count hydrogen atoms
        h_count = sum(1 for symbol in atoms_with_proton.get_chemical_symbols() if symbol == "H")
        assert h_count == 1

    # Check that locations are strings
    for location in locations:
        assert isinstance(location, str)


def test_add_two_protons():
    """Test adding two protons to the structure."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute two Si atoms with Al
    atoms = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Get oxygens and silicons around both substituted sites
    oxygens = []
    silicons = []
    for idx in CHA_DOUBLE_IDX:
        tmpo, tmps = get_os_and_ts(atoms, idx)
        oxygens.append(tmpo)
        silicons.append(tmps)

    # Add two protons
    traj, locations = add_two_protons(atoms, CHA_DOUBLE_IDX, oxygens, silicons, "CHA")

    # Check that structures and locations are returned
    assert isinstance(traj, list)
    assert isinstance(locations, list)
    assert len(traj) == len(locations)

    # Should get 16 structures (4x4 combinations)
    assert len(traj) == 16

    # Check that each structure has two protons added
    for atoms_with_protons in traj:
        assert "H" in atoms_with_protons.get_chemical_symbols()
        assert len(atoms_with_protons) == len(atoms) + 2
        # Count hydrogen atoms
        h_count = sum(1 for symbol in atoms_with_protons.get_chemical_symbols() if symbol == "H")
        assert h_count == 2

    # Check that locations are strings with proper format
    for location in locations:
        assert isinstance(location, str)
        assert "-" in location  # Should be in format "site1-site2"


def test_add_one_proton_with_path():
    """Test adding one proton with path parameter."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Get oxygens and silicons
    oxygens, silicons = get_os_and_ts(atoms, CHA_SINGLE_IDX)

    # Add one proton with path (but don't actually save files)
    traj, locations = add_one_proton(atoms, CHA_SINGLE_IDX, oxygens, silicons, "CHA", path=None)

    # Basic checks
    assert len(traj) == 4
    assert len(locations) == 4


def test_add_two_protons_with_path():
    """Test adding two protons with path parameter."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute two Si atoms with Al
    atoms = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Get oxygens and silicons
    oxygens = []
    silicons = []
    for idx in CHA_DOUBLE_IDX:
        tmpo, tmps = get_os_and_ts(atoms, idx)
        oxygens.append(tmpo)
        silicons.append(tmps)

    # Add two protons with path (but don't actually save files)
    traj, locations = add_two_protons(atoms, CHA_DOUBLE_IDX, oxygens, silicons, "CHA", path=None)

    # Basic checks
    assert len(traj) == 16
    assert len(locations) == 16
