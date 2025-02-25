import logging
import os.path
import sys

from .helpers import SeqType
from .modes import TrimmingMode
from .settings import DEFAULT_AA_GAP_CHARS

logger = logging.getLogger(__name__)


def process_args(args) -> dict:
    """
    Process args from argparser and set defaults
    """
    input_file = args.input
    output_file = args.output or f"{input_file}.clipkit"

    if not os.path.isfile(input_file):
        logger.warning("Input file does not exist")
        sys.exit()

    if input_file == output_file:
        logger.warning("Input and output files can't have the same name.")
        sys.exit()

    # assign optional arguments
    complement = args.complementary or False
    codon = args.codon or False
    mode = TrimmingMode(args.mode) if args.mode else TrimmingMode.smart_gap
    gaps = float(args.gaps) if args.gaps is not None else 0.9
    gap_characters = (
        [c for c in args.gap_characters] if args.gap_characters is not None else None
    )
    auxiliary_file = args.auxiliary_file
    use_log = args.log or False
    quiet = args.quiet or False
    sequence_type = SeqType(args.sequence_type.lower()) if args.sequence_type else None

    if codon and mode == TrimmingMode.c3:
        logger.warning(
            "C3 and codon-based trimming are incompatible.\nCodon-based trimming removes whole codons while C3 removes every third codon position."
        )
        sys.exit()

    ends_only = args.ends_only or False

    return dict(
        input_file=input_file,
        output_file=output_file,
        input_file_format=args.input_file_format,
        output_file_format=args.output_file_format,
        auxiliary_file=auxiliary_file,
        codon=codon,
        sequence_type=sequence_type,
        complement=complement,
        gaps=gaps,
        gap_characters=gap_characters,
        mode=mode,
        use_log=use_log,
        quiet=quiet,
        ends_only=ends_only,
    )
