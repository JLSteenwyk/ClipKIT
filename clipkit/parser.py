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
from .version import __version__

def create_parser():
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            f"""\
              _____ _ _       _  _______ _______ 
             / ____| (_)     | |/ /_   _|__   __|
            | |    | |_ _ __ | ' /  | |    | |   
            | |    | | | '_ \|  <   | |    | |   
            | |____| | | |_) | . \ _| |_   | |   
             \_____|_|_| .__/|_|\_\_____|  |_|   
                       | |                       
                       |_|  

        Version: {__version__}
        Citation: Steenwyk et al. 2020, PLOS Biology. doi: 10.1371/journal.pbio.3001007
        https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007

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

        -m, --modes <smart-gap,                     trimming mode 
                    gappy,                          (default: smart-gap)
                    kpic,
                    kpic-smart-gap,           
                    kpic-gappy,                
                    kpi,
                    kpi-smart-gap,
                    kpi-gappy>                      
                                                    
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
            smart-gap: dynamic determination of gaps threshold
            gappy: trim sites that are greater than the gaps threshold
            kpic: keeps parismony informative and constant sites
            kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
            kpic-gappy: a combination of kpic- and gappy-based trimming
            kpi: keep only parsimony informative sites
            kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
            kpi-gappy: a combination of kpi- and gappy-based trimming

        Gaps
            Positions with gappyness greater than threshold will be trimmed. 
            Must be between 0 and 1. (Default: 0.9). This argument is ignored
            when using the kpi and kpic mdoes of trimming as well as an 
            iteration of trimming that uses smart-gap.

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
        "-v", "--version", action="version", version=f"clipkit {__version__}", help=SUPPRESS,
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
