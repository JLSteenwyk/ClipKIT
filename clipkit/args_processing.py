
import os.path
import logging
import sys

from .modes import TrimmingMode

logger = logging.getLogger(__name__)

####################################################################
### Function to process arguments                                ###
### This function processes arguments                            ###
### subfunctions to trim the input file  						 ###
####################################################################
def process_args(args):
    inFile = args.input
    outFile = args.output or f"{inFile}.clipkit"

    # check that input file exists
    if not os.path.isfile(inFile):
        logger.info("Input file does not exist")
        sys.exit()

    if inFile == outFile:
        logger.info("Input and output files can't have the same name.")
        sys.exit()

    ## assign optional arguments
    mode = args.mode or TrimmingMode.gappy
    use_log = args.log or False
    gaps = float(args.gaps) if args.gaps else 0.9
    complement = args.complementary or False

    inFileFormat = args.input_file_format
    outFileFormat = args.output_file_format

    if use_log:
        # write INFO level logging to file for user
        logger.setLevel(logging.INFO)
        log_file_name = f"{outFile}.log"
        fh = logging.FileHandler(log_file_name, mode="w")
        fh.setLevel(logging.INFO)
        logger.addHandler(fh)

    return (
        inFile,
        outFile,
        inFileFormat,
        outFileFormat,
        gaps,
        complement,
        mode
    )

####################################################################
### END process_args Function    					             ###
####################################################################