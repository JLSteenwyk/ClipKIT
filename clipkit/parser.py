import sys
import textwrap

from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .modes import TrimmingMode
from .files import FileFormat

def create_parser():
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            """\
              _____ _ _       _  _______ _______ 
             / ____| (_)     | |/ /_   _|__   __|
            | |    | |_ _ __ | ' /  | |    | |   
            | |    | | | '_ \|  <   | |    | |   
            | |____| | | |_) | . \ _| |_   | |   
             \_____|_|_| .__/|_|\_\_____|  |_|   
                       | |                       
                       |_|  

        Citation: Steenwyk et al. bioRxiv.
        https://www.biorxiv.org/content/10.1101/2020.06.08.140384v1

        ClipKIT trims multiple sequence alignments and maintains phylogenetically informative sites.

        Usage: clipkit <input> [optional arguments]
        """
        ),
    )

    # if no arguments are given, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    ## required arguments
    required = parser.add_argument_group(
        "required arguments",
        description=textwrap.dedent(
            """\
        <input>                                     input file
                                                    (must be the first argument)
        """
        ),
    )

    required.add_argument("input", type=str, help=SUPPRESS)

    ## optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -o, --output <output_file_name>             output file name 
                                                    (default: input file named with '.clipkit' suffix)

        -m, --modes <gappy,                         trimming mode 
                    kpic (alias: medium),           (default: gappy)
                    kpic-gappy (alias: medium-gappy),                
                    kpi (alias: heavy),
                    kpi-gappy (alias: heavy-gappy)>                      
                                                    
        -g, --gaps <threshold of gaps>              specifies gaps threshold
                                                    (default: 0.9)

        -if, --input_file_format <file_format>      specifies input file format
                                                    (default: auto-detect)    

        -of, --output_file_format <file_format>     specifies output file format
                                                    (default: same as input file format)

        -l, --log                                   creates a log file
                                                    (input file named with '.log' suffix)

        -c, --complementary                         creates complementary alignment of trimmed sequences
                                                    (input file named with '.log' suffix)
 
        -h, --help                                  help message
        -v, --version                               print version

        -------------------------------------
        | Detailed explanation of arguments | 
        -------------------------------------
        Modes
            gappy: trim sites that are greater than the gaps threshold
            kpic (alias: medium): keeps parismony informative and constant sites
            kpic-gappy (alias: medium-gappy): a combination of kpic- and gappy-based trimming
            kpi (alias: heavy): keep only parsimony informative sites
            kpi-gappy (alias: heavy-gappy): a combination of kpi- and gappy-based trimming

        Gaps
            Positions with gappyness greater than threshold will be trimmed. 
            Must be between 0 and 1. (Default: 0.9). This argument is ignored
            when using the kpi mode of trimming.

        Input and output file formats
            Supported input and output files include:
            fasta, clustal, maf, mauve, phylip, phylip-sequential, 
            phylip-relaxed, and stockholm

        Log
            Creates a log file that summarizes the characteristics of each position.
            The log file has four columns.
            - Column 1 is the position in the alignment (starting at 1), 
            - Column 2 reports if the site was trimmed or kept (trim and keep, respectively),
            - Column 3 reports if the site is a parsimony informative site or not (PI and nPI, respectively), or
              a constant site or not (Const, nConst), or neither (nConst, nPI)
            - Column 4 reports the gappyness of the the position (number of gaps / entries in alignment)

        Complementary
            Creates an alignment file of only the trimmed sequences
        """
        ),
    )

    optional.add_argument("-o", "--output", help=SUPPRESS, metavar="output")

    mode_choices = [mode.value for mode in TrimmingMode]
    optional.add_argument(
        "-m", "--mode", help=SUPPRESS, nargs="?", choices=mode_choices,
    )

    optional.add_argument(
        "-h", "--help", action="help", help=SUPPRESS,
    )

    optional.add_argument(
        "-v", "--version", action="version", version="%(prog)s 0.1", help=SUPPRESS,
    )

    optional.add_argument(
        "-g",
        "--gaps",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="fraction of gaps",
    )

    file_format_choices = [file_format.value.lower() for file_format in FileFormat]
    optional.add_argument(
        "-if",
        "--input_file_format",
        type=str,
        required=False,
        choices=file_format_choices,
        help=SUPPRESS,
        metavar="",
    )

    optional.add_argument(
        "-of",
        "--output_file_format",
        type=str,
        required=False,
        choices=file_format_choices,
        help=SUPPRESS,
        metavar="",
    )

    optional.add_argument(
        "-l", "--log", action="store_true", required=False, help=SUPPRESS,
    )

    optional.add_argument(
        "-c", "--complementary", action="store_true", required=False, help=SUPPRESS,
    )

    return parser
