from __future__ import annotations

import pkg_resources
from ase.db import connect

path = "frameworks.db"
filepath = pkg_resources.resource_filename(__name__, path)


def get_framework(code):
    db = connect(filepath)
    return db.get_atoms(fw=code)


def get_ring_sizes(code):
    db = connect(filepath)
    return db.get(fw=code).data.rings


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
    return [row.fw for row in db.select()]
