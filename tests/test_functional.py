from pathlib import Path

from ase.io import read, write

from zse.substitute import exchange_unique_T_sites
from zse.utilities import make_iza_zeolite

REF_DATA = Path(__file__, "data")


def test_all_silica(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    atoms = make_iza_zeolite("CHA")
    write("CHA.cif", atoms)

    assert read("CHA.cif") == read(f"{REF_DATA}/CHA.cif")


def test_exchange_unique(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    zeolite = read(f"{REF_DATA}/CHA.cif")
    exchanged_zeolites = exchange_unique_T_sites(zeolite, "CHA", "B", "Na")
    assert len(exchanged_zeolites) == 4
    write("CHA_B_Na_3.cif", exchanged_zeolites[3])
    assert read("CHA_B_Na_3.cif") == read(f"{REF_DATA}/MOR_B_Na_3.cif")
