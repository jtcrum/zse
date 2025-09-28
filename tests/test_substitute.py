from zse.collections import framework
from zse.substitute import nest, tsub

CHA_SINGLE_IDX = 101
CHA_DOUBLE_IDX = [98, 101]


def test_tsub_single_substitution():
    """Test substituting a single atom."""
    # Create a CHA framework
    atoms = framework("CHA")
    original_symbols = atoms.get_chemical_symbols()

    # Substitute one Si atom with Al
    atoms_sub = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Check that atoms object is returned
    assert hasattr(atoms_sub, "get_chemical_symbols")

    # Check that the number of atoms hasn't changed
    assert len(atoms_sub) == len(atoms)

    # Check that the substitution occurred
    new_symbols = atoms_sub.get_chemical_symbols()
    assert new_symbols[CHA_SINGLE_IDX] == "Al"
    assert original_symbols[CHA_SINGLE_IDX] == "Si"  # Should have been Si originally

    # Check that only one atom was changed
    differences = sum(
        1 for old, new in zip(original_symbols, new_symbols, strict=True) if old != new
    )
    assert differences == 1

    # Check that original atoms object is not modified
    assert atoms.get_chemical_symbols()[CHA_SINGLE_IDX] == "Si"


def test_tsub_multiple_substitution():
    """Test substituting multiple atoms."""
    # Create a CHA framework
    atoms = framework("CHA")
    original_symbols = atoms.get_chemical_symbols()

    # Substitute two Si atoms with Al
    atoms_sub = tsub(atoms, CHA_DOUBLE_IDX, "Al")

    # Check that the number of atoms hasn't changed
    assert len(atoms_sub) == len(atoms)

    # Check that both substitutions occurred
    new_symbols = atoms_sub.get_chemical_symbols()
    for idx in CHA_DOUBLE_IDX:
        assert new_symbols[idx] == "Al"
        assert original_symbols[idx] == "Si"  # Should have been Si originally

    # Check that exactly two atoms were changed
    differences = sum(
        1 for old, new in zip(original_symbols, new_symbols, strict=True) if old != new
    )
    assert differences == 2


def test_tsub_different_elements():
    """Test substituting with different elements."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test substitution with different elements
    test_elements = ["Al", "B", "Ga", "Fe"]

    for element in test_elements:
        atoms_sub = tsub(atoms, CHA_SINGLE_IDX, element)
        new_symbols = atoms_sub.get_chemical_symbols()
        assert new_symbols[CHA_SINGLE_IDX] == element


def test_tsub_list_input():
    """Test tsub with list input for index."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test with single-element list
    atoms_sub1 = tsub(atoms, [CHA_SINGLE_IDX], "Al")
    assert atoms_sub1.get_chemical_symbols()[CHA_SINGLE_IDX] == "Al"

    # Test with multi-element list
    atoms_sub2 = tsub(atoms, CHA_DOUBLE_IDX, "Al")
    new_symbols = atoms_sub2.get_chemical_symbols()
    for idx in CHA_DOUBLE_IDX:
        assert new_symbols[idx] == "Al"


def test_tsub_preserves_structure():
    """Test that tsub preserves the structure (positions, cell, etc.)."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Get original properties
    original_positions = atoms.get_positions()
    original_cell = atoms.cell
    original_pbc = atoms.pbc

    # Substitute one atom
    atoms_sub = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Check that structure properties are preserved
    new_positions = atoms_sub.get_positions()
    assert len(new_positions) == len(original_positions)

    # Positions should be identical
    assert abs(new_positions - original_positions).max() < 1e-10

    # Cell should be identical
    assert abs(atoms_sub.cell.array - original_cell.array).max() < 1e-10

    # PBC should be identical
    assert all(atoms_sub.pbc == original_pbc)


def test_nest_basic():
    """Test basic nest defect creation."""
    # Create a CHA framework
    atoms = framework("CHA")
    original_length = len(atoms)

    # Create nest defect
    nest_traj = nest(atoms, CHA_SINGLE_IDX)

    # Check that a list of structures is returned
    assert isinstance(nest_traj, list)

    # Should get 6 structures (different proton configurations)
    assert len(nest_traj) == 6

    # Check each structure
    for nest_atoms in nest_traj:
        # Should have hydrogen atoms added and one Si removed
        symbols = nest_atoms.get_chemical_symbols()
        assert "H" in symbols

        # Count hydrogen atoms (should be 4)
        h_count = sum(1 for symbol in symbols if symbol == "H")
        assert h_count == 4

        # Total atoms should be original - 1 (removed Si) + 4 (added H) = original + 3
        assert len(nest_atoms) == original_length + 3

        # Check that the original T-site atom is removed
        # The nest structure should not have the original atom at that index
        # (structure is modified so direct index comparison may not work,
        # but we can check that one Si atom is missing)
        original_si_count = sum(1 for symbol in atoms.get_chemical_symbols() if symbol == "Si")
        new_si_count = sum(1 for symbol in symbols if symbol == "Si")
        assert new_si_count == original_si_count - 1


def test_nest_preserves_cell():
    """Test that nest defect preserves the cell parameters."""
    # Create a CHA framework
    atoms = framework("CHA")
    original_cell = atoms.cell

    # Create nest defect
    nest_traj = nest(atoms, CHA_SINGLE_IDX)

    # Check that cell is preserved in all structures
    for nest_atoms in nest_traj:
        assert abs(nest_atoms.cell.array - original_cell.array).max() < 1e-10


def test_nest_different_indices():
    """Test nest creation at different T-site indices."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test different T-site indices
    test_indices = [50, 75, 100]

    for index in test_indices:
        # Make sure the index corresponds to a Si atom
        if atoms[index].symbol == "Si":
            nest_traj = nest(atoms, index)

            # Should always get 6 structures
            assert len(nest_traj) == 6

            # Each should have the same properties
            for nest_atoms in nest_traj:
                symbols = nest_atoms.get_chemical_symbols()
                h_count = sum(1 for symbol in symbols if symbol == "H")
                assert h_count == 4


def test_tsub_with_oxygen():
    """Test that tsub can work with oxygen atoms too."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Find an oxygen atom
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "O"]
    assert len(oxygen_indices) > 0

    # Substitute an oxygen with fluorine
    oxy_index = oxygen_indices[0]
    atoms_sub = tsub(atoms, oxy_index, "F")

    # Check substitution
    assert atoms_sub.get_chemical_symbols()[oxy_index] == "F"


def test_tsub_returns_copy():
    """Test that tsub returns a copy and doesn't modify original."""
    # Create a CHA framework
    atoms = framework("CHA")
    original_symbol = atoms[CHA_SINGLE_IDX].symbol

    # Substitute atom
    atoms_sub = tsub(atoms, CHA_SINGLE_IDX, "Al")

    # Check that original is unchanged
    assert atoms[CHA_SINGLE_IDX].symbol == original_symbol

    # Check that copy is changed
    assert atoms_sub[CHA_SINGLE_IDX].symbol == "Al"

    # Check that they are different objects
    assert atoms is not atoms_sub
