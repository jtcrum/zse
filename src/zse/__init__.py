from __future__ import annotations

from typing import Any

from ase.atoms import Atoms


def atoms_as_dict(s: Atoms) -> dict[str, Any]:
    from ase.io.jsonio import encode

    # Uses Monty's MSONable spec
    # Normally, we would want to this to be a wrapper around atoms.todict() with @module and
    # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() does not currently
    # work properly with constraints.
    return {"@module": "ase.atoms", "@class": "Atoms", "atoms_json": encode(s)}


def atoms_from_dict(d: dict[str, Any]) -> Atoms:
    from ase.io.jsonio import decode

    # Uses Monty's MSONable spec
    # Normally, we would want to have this be a wrapper around atoms.fromdict()
    # that just ignores the @module/@class key-value pairs. However, atoms.todict()/atoms.fromdict()
    # does not currently work properly with constraints.
    return decode(d["atoms_json"])


# Make Atoms MSONable
Atoms.as_dict = atoms_as_dict
Atoms.from_dict = atoms_from_dict
