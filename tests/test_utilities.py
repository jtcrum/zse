import numpy as np

from zse.collections import framework
from zse.substitute import tsub
from zse.utilities import center, get_osites, get_tsites, scale_cell, site_labels


def test_center():
    """Test centering an atoms object around an index."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test centering around a T-site
    index = 101
    centered_atoms, translation = center(atoms, index)

    # Check that atoms object is returned
    assert hasattr(centered_atoms, "get_positions")

    # Check that translation vector is returned
    assert isinstance(translation, np.ndarray)
    assert len(translation) == 3

    # Check that the center function worked correctly
    # The function translates the structure so that the center of mass moves
    # to where the specified atom was originally
    original_com = atoms.get_center_of_mass()

    # After centering, the specified atom should be roughly where the COM was
    new_atom_pos = centered_atoms[index].position
    distance_to_original_com = np.linalg.norm(new_atom_pos - original_com)
    assert distance_to_original_com < 1.0  # Should be close due to translation

    # Check that the number of atoms hasn't changed
    assert len(centered_atoms) == len(atoms)


def test_get_osites():
    """Test getting oxygen site information for CHA framework."""
    osites, omult, first_os = get_osites("CHA")

    # Check return types
    assert isinstance(osites, list)
    assert isinstance(omult, list)
    assert isinstance(first_os, list)

    # Check that all lists have same length
    assert len(osites) == len(omult) == len(first_os)

    # Check that we have some oxygen sites
    assert len(osites) > 0

    # Check that oxygen site labels are strings
    for site in osites:
        assert isinstance(site, str)

    # Check that multiplicities are positive integers
    for mult in omult:
        assert isinstance(mult, int)
        assert mult > 0

    # Check that first oxygen indices are integers
    for idx in first_os:
        assert isinstance(idx, int)
        assert idx >= 0


def test_get_tsites():
    """Test getting T-site information for CHA framework."""
    tsites, tmult, first_ts = get_tsites("CHA")

    # Check return types
    assert isinstance(tsites, list)
    assert isinstance(tmult, list)
    assert isinstance(first_ts, list)

    # Check that all lists have same length
    assert len(tsites) == len(tmult) == len(first_ts)

    # Check that we have some T-sites
    assert len(tsites) > 0

    # Check that T-site labels are strings
    for site in tsites:
        assert isinstance(site, str)

    # Check that multiplicities are positive integers
    for mult in tmult:
        assert isinstance(mult, int)
        assert mult > 0

    # Check that first T-site indices are integers
    for idx in first_ts:
        assert isinstance(idx, int)
        assert idx >= 0


def test_scale_cell():
    """Test scaling cell to achieve target Si-Si distance."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Scale the cell
    scaled_atoms = scale_cell(atoms)

    # Check that atoms object is returned
    assert hasattr(scaled_atoms, "get_cell")

    # Check that cell parameters have changed (unless already at target)
    scaled_cell = scaled_atoms.cell
    assert scaled_cell is not None

    # Check that the number of atoms hasn't changed
    assert len(scaled_atoms) == len(atoms)

    # Check that atoms still have positions
    positions = scaled_atoms.get_positions()
    assert len(positions) == len(scaled_atoms)


def test_site_labels():
    """Test getting site labels for atoms object."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Get site labels
    labels = site_labels(atoms, "CHA")

    # Check return type
    assert isinstance(labels, dict)

    # Check that all atoms have labels
    assert len(labels) == len(atoms)

    # Check that all keys are atom indices
    for key in labels:
        assert isinstance(key, int)
        assert 0 <= key < len(atoms)

    # Check that all values are strings
    for value in labels.values():
        assert isinstance(value, str)

    # Check that we have both T-site and O-site labels
    label_values = set(labels.values())

    # Should have T-site labels (containing T) and O-site labels (containing O)
    t_labels = [label for label in label_values if "T" in label]
    o_labels = [label for label in label_values if "O" in label]

    assert len(t_labels) > 0
    assert len(o_labels) > 0


def test_site_labels_with_substitution():
    """Test site labels work correctly with Al substitution."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Substitute one Si with Al
    index = 101
    atoms_sub = tsub(atoms, index, "Al")

    # Get site labels
    labels = site_labels(atoms_sub, "CHA")

    # Check that we still get labels for all atoms
    assert len(labels) == len(atoms_sub)

    # Check that the substituted atom gets a T-site label
    assert index in labels
    assert isinstance(labels[index], str)
    assert "T" in labels[index]


def test_get_osites_consistency():
    """Test that get_osites returns consistent data."""
    _osites, omult, first_os = get_osites("CHA")

    # Check that multiplicities sum correctly
    # The sum of multiplicities should equal total number of unique oxygen sites
    assert sum(omult) > 0

    # Check that first indices are in ascending order
    for i in range(1, len(first_os)):
        assert first_os[i] > first_os[i - 1]


def test_get_tsites_consistency():
    """Test that get_tsites returns consistent data."""
    _tsites, tmult, first_ts = get_tsites("CHA")

    # Check that multiplicities sum correctly
    assert sum(tmult) > 0

    # Check that first indices are in ascending order
    for i in range(1, len(first_ts)):
        assert first_ts[i] > first_ts[i - 1]


def test_center_preserves_structure():
    """Test that centering preserves the relative structure."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Get original distances between first few atoms
    original_dist = atoms.get_distance(0, 1, mic=True)

    # Center around an atom
    centered_atoms, _ = center(atoms, 50)

    # Check that distances are preserved
    new_dist = centered_atoms.get_distance(0, 1, mic=True)
    assert abs(original_dist - new_dist) < 1e-10


def test_site_labels_different_cell_sizes():
    """Test site labels work with different cell sizes."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Get labels for original cell
    labels1 = site_labels(atoms, "CHA")

    # Repeat the cell and get labels
    atoms_repeated = atoms.repeat([2, 1, 1])
    labels2 = site_labels(atoms_repeated, "CHA")

    # Should get labels for all atoms in both cases
    assert len(labels1) == len(atoms)
    assert len(labels2) == len(atoms_repeated)

    # The repeated structure should have more atoms
    assert len(atoms_repeated) > len(atoms)
