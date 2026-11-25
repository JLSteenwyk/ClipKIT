from enum import Enum


class FileFormat(Enum):
    # Preserve the original import paths for serialized enum values.
    __module__ = "clipkit.files"

    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_sequential = "phylip_sequential"
    phylip_relaxed = "phylip_relaxed"
    stockholm = "stockholm"
    ecomp = "ecomp"


class SeqType(Enum):
    __module__ = "clipkit.helpers"

    aa = "aa"
    nt = "nt"
