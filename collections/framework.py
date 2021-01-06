from ase.db import connect
import pkg_resources

path = 'frameworks.db'
filepath = pkg_resources.resource_filename(__name__,path)

def framework(code):
    db = connect(filepath)
    atoms = db.get_atoms(fw = code)
    return atoms

def get_fwrings(code):
    db = connect(filepath)
    fw_rings = db.get(fw = code).data.rings
    return fw_rings

def get_tsites(code):
    db = connect(filepath)
    tsites = db.get(fw=code).data.tsites
    tmult = db.get(fw=code).data.tmult
    return tsites,tmult

def get_osites(code):
    db = connect(fielpath)
    osites = db.get(fw=code).data.osites
    omult = db.get(fw=code).data.omult
    return osites,omult

def get_all_fws():
    db = connect(filepath)
    fws = []
    for row in db.select():
        fws.append(row.fw)
    return fws
