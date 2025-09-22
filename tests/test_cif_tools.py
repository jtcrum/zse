import os

from ase import Atoms

from zse.cif_tools import cif_site_labels, read_cif

# Path to test CIF file
TEST_CIF_PATH = os.path.join(os.path.dirname(__file__), "cif_example", "BEC.cif")


def test_read_cif():
    """Test reading a CIF file and extracting site information."""
    # Read the test CIF file
    atoms, ts, tm, tinds, os_sites, om, oinds = read_cif(TEST_CIF_PATH)

    # Check that an Atoms object is returned
    assert isinstance(atoms, Atoms)

    # Check that T-site information is returned
    assert isinstance(ts, list)
    assert isinstance(tm, list)
    assert isinstance(tinds, list)
    assert len(ts) == len(tm) == len(tinds)

    # Check that O-site information is returned
    assert isinstance(os_sites, list)
    assert isinstance(om, list)
    assert isinstance(oinds, list)
    assert len(os_sites) == len(om) == len(oinds)

    # Check that we have some T-sites and O-sites
    assert len(ts) > 0
    assert len(os_sites) > 0

    # Check that T-site indices are valid
    for idx in tinds:
        assert 0 <= idx < len(atoms)

    # Check that O-site indices are valid
    for idx in oinds:
        assert 0 <= idx < len(atoms)

    # Check that T-site labels are strings
    for label in ts:
        assert isinstance(label, str)

    # Check that O-site labels are strings
    for label in os_sites:
        assert isinstance(label, str)

    # Check that multiplicities are positive integers
    for mult in tm:
        assert isinstance(mult, int)
        assert mult > 0

    for mult in om:
        assert isinstance(mult, int)
        assert mult > 0


def test_cif_site_labels():
    """Test generation of site label mapping from CIF file."""
    # Get site labels mapping
    labels = cif_site_labels(TEST_CIF_PATH)

    # Check that a dictionary is returned
    assert isinstance(labels, dict)

    # Check that the dictionary is not empty
    assert len(labels) > 0

    # Check that all keys are integers (atom indices)
    for key in labels:
        assert isinstance(key, int)

    # Check that all values are strings (site labels)
    for value in labels.values():
        assert isinstance(value, str)

    # Read atoms to check index validity
    atoms, _ts, _tm, _tinds, _os, _om, _oinds = read_cif(TEST_CIF_PATH)

    # Check that all indices in labels are valid atom indices
    for idx in labels:
        assert 0 <= idx < len(atoms)

    # Check that we have labels for both T-sites and O-sites
    label_values = set(labels.values())
    # Should have some T-site labels (containing T) and O-site labels (containing O)
    t_labels = [label for label in label_values if "T" in label]
    o_labels = [label for label in label_values if "O" in label]

    assert len(t_labels) > 0
    assert len(o_labels) > 0


def test_cif_file_exists():
    """Test that the test CIF file exists."""
    assert os.path.exists(TEST_CIF_PATH), f"Test CIF file not found at {TEST_CIF_PATH}"


def test_read_cif_atoms_properties():
    """Test that the atoms object from CIF has expected properties."""
    atoms, _ts, _tm, _tinds, _os, _om, _oinds = read_cif(TEST_CIF_PATH)

    # Check that atoms object has expected properties
    assert len(atoms) > 0
    assert atoms.cell is not None
    assert len(atoms.cell) == 3

    # Check that we have both Si and O atoms
    symbols = atoms.get_chemical_symbols()
    assert "Si" in symbols or "T" in symbols  # T-sites might be labeled as T
    assert "O" in symbols

    # Check that positions are available
    positions = atoms.get_positions()
    assert len(positions) == len(atoms)
    assert positions.shape == (len(atoms), 3)


def test_cif_site_consistency():
    """Test consistency between read_cif and cif_site_labels."""
    # Read CIF file
    _atoms, ts, tm, tinds, os_sites, om, oinds = read_cif(TEST_CIF_PATH)

    # Get site labels mapping
    labels = cif_site_labels(TEST_CIF_PATH)

    # Check that T-site labels are consistent
    for i, t_label in enumerate(ts):
        for j in range(tm[i]):
            atom_idx = tinds[i] + j
            assert atom_idx in labels
            assert labels[atom_idx] == t_label

    # Check that O-site labels are consistent
    for i, o_label in enumerate(os_sites):
        for j in range(om[i]):
            atom_idx = oinds[i] + j
            assert atom_idx in labels
            assert labels[atom_idx] == o_label
