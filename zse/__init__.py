"""Zeolite Simulation Environment (ZSE).

A Python package for automating zeolite structure generation and characterization
for computational chemistry.
"""

__version__ = "0.0.1"
__author__ = "Jerry T. Crum, Justin R. Crum"

# Import main modules
from . import cation, protonate, rings, substitute, tpairs, utilities
from .collections import framework

__all__ = [
    "cation",
    "framework",
    "protonate",
    "rings",
    "substitute",
    "tpairs",
    "utilities",
]