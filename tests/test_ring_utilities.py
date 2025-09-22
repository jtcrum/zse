import networkx as nx

from zse.collections import framework
from zse.ring_utilities import atoms_to_graph, remove_dups, vertex_order


def test_atoms_to_graph():
    """Test converting atoms to graph representation."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test atoms_to_graph function
    index = 101
    max_ring = 12
    graph, large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)

    # Check return types
    assert isinstance(graph, nx.Graph)
    assert hasattr(large_atoms, "get_positions")
    assert isinstance(repeat, list)

    # Check that graph has nodes and edges
    assert len(graph.nodes) > 0
    assert len(graph.edges) > 0

    # Check that repeat is a valid list
    assert len(repeat) == 3
    for r in repeat:
        assert isinstance(r, int)
        assert r >= 1

    # Check that large_atoms is larger than or equal to original
    assert len(large_atoms) >= len(atoms)


def test_remove_dups():
    """Test removing duplicate paths."""
    # Create some test paths with duplicates
    paths = [
        [1, 2, 3, 4],
        [4, 3, 2, 1],  # Reverse of first
        [1, 2, 3, 4],  # Exact duplicate
        [5, 6, 7, 8],  # Different path
        [8, 7, 6, 5],  # Reverse of different path
    ]

    # Remove duplicates
    unique_paths = remove_dups(paths)

    # Should have fewer paths than input
    assert len(unique_paths) <= len(paths)

    # Should have at least some paths
    assert len(unique_paths) > 0

    # Check that returned value is a list
    assert isinstance(unique_paths, list)


def test_vertex_order():
    """Test vertex ordering functionality."""
    # Create a CHA framework and get a ring
    atoms = framework("CHA")

    # Import here to avoid circular imports
    from zse.rings import get_rings  # noqa: PLC0415

    # Get rings for testing
    _ring_list, paths, _ring_atoms, atoms = get_rings(atoms, 101, validation="crum")

    if len(paths) > 0:
        # Test with the first ring
        test_paths = paths[:2]
        print(test_paths)

        # Test vertex_order function
        ordered, _new_rings = vertex_order(test_paths)

        # Check that result is a string
        assert isinstance(ordered, str)

        # Should contain ring size information
        assert any(str(size) in ordered for size in [4, 6, 8, 10, 12])


def test_atoms_to_graph_different_indices():
    """Test atoms_to_graph with different indices."""
    # Create a CHA framework
    atoms = framework("CHA")

    # Test with different T-site indices
    test_indices = [50, 75, 100]
    max_ring = 8

    for index in test_indices:
        # Make sure the index corresponds to a Si atom
        if atoms[index].symbol == "Si":
            graph, _large_atoms, repeat = atoms_to_graph(atoms, index, max_ring)

            # Basic checks
            assert isinstance(graph, nx.Graph)
            assert len(graph.nodes) > 0
            assert len(graph.edges) > 0
            assert len(repeat) == 3


def test_remove_dups_empty():
    """Test remove_dups with empty input."""
    # Test with empty list
    result = remove_dups([])
    assert result == []


def test_remove_dups_single():
    """Test remove_dups with single path."""
    # Test with single path
    paths = [[1, 2, 3, 4]]
    result = remove_dups(paths)
    assert len(result) == 1
    assert result[0] == [1, 2, 3, 4]


def test_atoms_to_graph_properties():
    """Test that atoms_to_graph creates a proper graph."""
    # Create a CHA framework
    atoms = framework("CHA")

    index = 101
    max_ring = 10
    graph, large_atoms, _repeat = atoms_to_graph(atoms, index, max_ring)

    # Check graph properties
    assert graph.number_of_nodes() > 0
    assert graph.number_of_edges() > 0

    # Check that large_atoms has expected properties
    assert len(large_atoms) > 0
    positions = large_atoms.get_positions()
    assert positions.shape[0] == len(large_atoms)
    assert positions.shape[1] == 3
