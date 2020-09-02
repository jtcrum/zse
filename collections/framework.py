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
