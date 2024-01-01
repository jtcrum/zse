from pathlib import Path

from ase.io import read
from numpy.testing import assert_array_almost_equal, assert_array_equal

from zse.substitute import exchange_unique_T_sites
from zse.utilities import make_iza_zeolite

REF_DATA = Path(__file__).parent / Path("data")


def test_all_silica():
    atoms = make_iza_zeolite("CHA")
    ref_atoms = read(f"{REF_DATA}/CHA.cif")
    assert_array_equal(atoms.get_chemical_symbols(), ref_atoms.get_chemical_symbols())
    assert_array_almost_equal(atoms.get_positions(), ref_atoms.get_positions())
    assert_array_almost_equal(atoms.get_cell(), ref_atoms.get_cell())


def test_exchange_unique():
    exchanged_zeolites = exchange_unique_T_sites("CHA", "B", "Na")
    assert len(exchanged_zeolites) == 4

    atoms = exchanged_zeolites[3]
    ref_atoms = read(f"{REF_DATA}/MOR_B_Na_3.cif")
    assert_array_equal(atoms.get_chemical_symbols(), ref_atoms.get_chemical_symbols())
    assert_array_almost_equal(atoms.get_positions(), ref_atoms.get_positions())
    assert_array_almost_equal(atoms.get_cell(), ref_atoms.get_cell())
