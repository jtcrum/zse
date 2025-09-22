from zse.collections import framework
from zse.protonate import isolated, paired
from zse.substitute import tsub

CHA_SINGLE_IDX = 101
CHA_DOUBLE_IDX = [98, 101]


def test_isolated():
    """Test isolated protonation functionality."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Test isolated protonation
    traj, locations = isolated(atoms, CHA_SINGLE_IDX, "CHA")

    # Check that structures and locations are returned
    assert isinstance(traj, list)
    assert isinstance(locations, list)
    assert len(traj) == len(locations)
    assert len(traj) > 0

    # Check that each structure has a proton added
    for atoms_with_proton in traj:
        assert "H" in atoms_with_proton.get_chemical_symbols()
        assert len(atoms_with_proton) == len(atoms) + 1


def test_paired():
    """Test paired protonation functionality."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute two Si atoms with Al
    atoms = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Test paired protonation
    traj, locations = paired(atoms, CHA_DOUBLE_IDX, "CHA")

    # Check that structures and locations are returned
    assert isinstance(traj, list)
    assert isinstance(locations, list)
    assert len(traj) == len(locations)
    assert len(traj) > 0

    # Check that each structure has two protons added
    for atoms_with_protons in traj:
        assert "H" in atoms_with_protons.get_chemical_symbols()
        assert len(atoms_with_protons) == len(atoms) + 2
        # Count hydrogen atoms
        h_count = sum(1 for symbol in atoms_with_protons.get_chemical_symbols() if symbol == "H")
        assert h_count == 2


def test_isolated_with_path():
    """Test isolated protonation with path parameter."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si atom with Al
    atoms = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Test isolated protonation with path (but don't actually save files)
    traj, locations = isolated(atoms, CHA_SINGLE_IDX, "CHA", path=None)

    # Basic checks
    assert len(traj) > 0
    assert len(locations) > 0
    assert len(traj) == len(locations)


def test_paired_with_path():
    """Test paired protonation with path parameter."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute two Si atoms with Al
    atoms = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Test paired protonation with path (but don't actually save files)
    traj, locations = paired(atoms, CHA_DOUBLE_IDX, "CHA", path=None)

    # Basic checks
    assert len(traj) > 0
    assert len(locations) > 0
    assert len(traj) == len(locations)
