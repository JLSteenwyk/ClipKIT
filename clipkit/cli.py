import sys

from .exceptions import StopCodonValidationError
from .parser import create_parser


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    parser = create_parser()
    argv = sys.argv[1:] if argv is None else argv

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        return

    args = parser.parse_args(argv)

    from .args_processing import process_args
    from .clipkit import execute

    try:
        execute(**process_args(args))
    except StopCodonValidationError as exc:
        parser.error(str(exc))
