from ase.db import connect
import pkg_resources

path = 'frameworks.db'
filepath = pkg_resources.resource_filename(__name__,path)

def framework(code):
    db = connect(filepath)
    atoms = db.get_atoms(fw = code)
    return atoms
