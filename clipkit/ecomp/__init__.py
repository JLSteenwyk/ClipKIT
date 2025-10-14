"""Utilities for reading and writing Evolutionary Compression (.ecomp) archives."""

from .reader import read_ecomp
from .writer import write_ecomp

__all__ = ["read_ecomp", "write_ecomp"]
