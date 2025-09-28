import numpy as np
from ase import Atoms

from zse.collections.framework import (
    framework_db,
    get_all_fws_db,
    get_osites_db,
    get_ring_sizes_db,
    get_tsites_db,
)


def test_framework_db():
    """Test getting framework from database."""
    # Test with CHA framework
    atoms = framework_db("CHA")

    # Check that an Atoms object is returned
    assert isinstance(atoms, Atoms)

    # Check that it has the expected properties
    assert len(atoms) > 0
    assert atoms.cell is not None

    # Check that it contains Si and O atoms
    symbols = atoms.get_chemical_symbols()
    assert "Si" in symbols
    assert "O" in symbols


def test_get_ring_sizes_db():
    """Test getting ring sizes from database."""
    # Test with CHA framework
    ring_sizes = get_ring_sizes_db("CHA")

    # Check that a list is returned
    assert isinstance(ring_sizes, (np.ndarray, list))

    # Check that we get some ring sizes
    assert len(ring_sizes) > 0

    # Check that all ring sizes are integers
    for size in ring_sizes:
        assert isinstance(size, (np.int8, np.int16, np.int32, np.int64, int))
        assert size > 0


def test_get_tsites_db():
    """Test getting T-sites from database."""
    # Test with CHA framework
    tsites, tmult = get_tsites_db("CHA")

    # Check return types
    assert isinstance(tsites, list)
    assert isinstance(tmult, list)

    # Check that both lists have same length
    assert len(tsites) == len(tmult)

    # Check that we have some T-sites
    assert len(tsites) > 0

    # Check that T-site labels are strings
    for site in tsites:
        assert isinstance(site, str)

    # Check that multiplicities are positive integers
    for mult in tmult:
        assert isinstance(mult, int)
        assert mult > 0


def test_get_osites_db():
    """Test getting oxygen sites from database."""
    # Test with CHA framework
    osites, omult = get_osites_db("CHA")

    # Check return types
    assert isinstance(osites, list)
    assert isinstance(omult, list)

    # Check that both lists have same length
    assert len(osites) == len(omult)

    # Check that we have some oxygen sites
    assert len(osites) > 0

    # Check that oxygen site labels are strings
    for site in osites:
        assert isinstance(site, str)

    # Check that multiplicities are positive integers
    for mult in omult:
        assert isinstance(mult, int)
        assert mult > 0


def test_get_all_fws_db():
    """Test getting all framework codes from database."""
    # Get all framework codes
    fws = get_all_fws_db()

    # Check that a list is returned
    assert isinstance(fws, list)

    # Check that we have some frameworks
    assert len(fws) > 0

    # Check that all codes are strings
    for fw in fws:
        assert isinstance(fw, str)

    # Check that CHA is in the list (since we know it exists)
    assert "CHA" in fws


def test_framework_consistency():
    """Test consistency between different database functions."""
    # Get framework and its properties
    atoms = framework_db("CHA")
    _tsites, tmult = get_tsites_db("CHA")
    _osites, omult = get_osites_db("CHA")

    # Count actual atoms
    symbols = atoms.get_chemical_symbols()
    si_count = sum(1 for symbol in symbols if symbol == "Si")
    o_count = sum(1 for symbol in symbols if symbol == "O")

    # Check that T-site multiplicities sum to number of Si atoms
    total_t_atoms = sum(tmult)
    assert total_t_atoms == si_count

    # Check that O-site multiplicities sum to number of O atoms
    total_o_atoms = sum(omult)
    assert total_o_atoms == o_count


def test_multiple_frameworks():
    """Test that multiple frameworks can be loaded."""
    # Get all available frameworks
    fws = get_all_fws_db()

    # Test a few different frameworks (if available)
    test_frameworks = ["CHA"]
    for fw in test_frameworks:
        if fw in fws:
            atoms = framework_db(fw)
            assert isinstance(atoms, Atoms)
            assert len(atoms) > 0

            # Test that we can get sites for this framework
            tsites, _tmult = get_tsites_db(fw)
            osites, _omult = get_osites_db(fw)

            assert len(tsites) > 0
            assert len(osites) > 0
