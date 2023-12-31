from ase.io import read

from zse.substitute import exchange_unique_T_sites
from zse.utilities import make_iza_zeolite


def test_all_silica():
    atoms = make_iza_zeolite("CHA")
    assert atoms == read("data/CHA.cif")


def test_exchange_unique():
    zeolite = read(f"MOR.cif")
    exchanged_zeolites = exchange_unique_T_sites(zeolite, "CHA", "B", "Na")
    assert len(exchanged_zeolites) == 4
    assert exchanged_zeolites[0] == read("data/MOR_B_Na_1.cif")
