from typing import TYPE_CHECKING

import numpy as np
import pkg_resources
from ase.db import connect

from zse.utilities import site_labels

if TYPE_CHECKING:
    from ase.atoms import Atoms
path = "frameworks.db"
filepath = pkg_resources.resource_filename(__name__, path)


def make_iza_zeolite(code: str) -> Atoms:
    """
    Make an idealized zeolite from an IZA code, populate the atoms.info
    dictionary with the framework name, and add a labels array to the
    atoms object.
    """
    zeolite = framework(code)
    labels = site_labels(zeolite, code)
    zeolite.set_array("labels", np.array(list(labels.values())))
    zeolite.info["framework"] = code
    return zeolite


def framework(code):
    db = connect(filepath)
    atoms = db.get_atoms(fw=code)
    return atoms


def get_ring_sizes(code):
    db = connect(filepath)
    fw_rings = db.get(fw=code).data.rings
    return fw_rings


def get_tsites(code):
    db = connect(filepath)
    tsites = db.get(fw=code).data.tsites
    tmult = db.get(fw=code).data.tmult
    return tsites, tmult


def get_osites(code):
    db = connect(filepath)
    osites = db.get(fw=code).data.osites
    omult = db.get(fw=code).data.omult
    return osites, omult


def get_all_fws():
    db = connect(filepath)
    fws = []
    for row in db.select():
        fws.append(row.fw)
    return fws
