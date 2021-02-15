import logging
import os.path
import sys

from .modes import TrimmingMode

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
    mode = TrimmingMode(args.mode) if args.mode else TrimmingMode.smart_gap
    gaps = float(args.gaps) if args.gaps is not None else 0.9
    use_log = args.log or False

    return dict(
        input_file=input_file,
        output_file=output_file,
        input_file_format=args.input_file_format,
        output_file_format=args.output_file_format,
        complement=complement,
        gaps=gaps,
        mode=mode,
        use_log=use_log,
    )
