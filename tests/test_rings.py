from zse.collections import framework
from zse.rings import get_ordered_vertex, get_rings, get_unique_rings

CHA_SINGLE_IDX = 101


def test_get_rings():
    # Create a CHA framework
    atoms = framework("CHA")

    # Get the rings for T-site
    ring_list, paths, ring_atoms, atoms = get_rings(atoms, CHA_SINGLE_IDX, validation="crum")

    # Check that the correct number of rings are returned
    assert len(ring_list) == 7
    assert len(paths) == 7
    assert len(ring_atoms) == 7

    # Check that the ring sizes are correct
    assert ring_list.count(4) == 3
    assert ring_list.count(6) == 1
    assert ring_list.count(8) == 2


def test_get_unique_rings():
    # Create a CHA framework
    atoms = framework("CHA")

    # Get the unique rings for the framework
    ring_list, paths, ring_atoms, atoms = get_unique_rings(
        atoms, [CHA_SINGLE_IDX], validation="crum"
    )

    # Check that the correct number of unique rings are returned
    assert len(ring_list) == 5
    assert len(paths) == 5
    assert len(ring_atoms) == 5
    assert 4 in ring_list
    assert 6 in ring_list
    assert 8 in ring_list


def test_get_ordered_vertex():
    # Create a CHA framework
    atoms = framework("CHA")

    # Get the ordered vertex for T-site 48
    ordered_vertex, _paths, _ring_atoms, atoms = get_ordered_vertex(atoms, CHA_SINGLE_IDX)

    # Check that the ordered vertex is correct
    assert "8" in ordered_vertex
    assert "6" in ordered_vertex
    assert "4" in ordered_vertex
