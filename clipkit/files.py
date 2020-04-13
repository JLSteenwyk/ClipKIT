from Bio import AlignIO

## Function to automatically determine the format of the alignment file
## and read in the alignment. Returns alignment object and fileFormat
def automatic_file_type_determination(
    inFile
    ):
    """
    Automatically determines what type of input file was used
    and reads in the alignment file

    Arguments
    ---------
    argv: inFile
        input file specified with -i, --input
    """

    # save list of different file formats
    fileFormats = ['fasta', 'clustal', 'maf', 'mauve',
        'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm']

    # loop through file formats and attempt to read in file in that format
    for fileFormat in fileFormats:
        try:
            alignment = AlignIO.read(open(inFile), fileFormat)
            return alignment, fileFormat
            break
        # the following exceptions refer to skipping over errors 
        # associated with reading the wrong input file
        except ValueError:
            continue
        except AssertionError:
            continue

## Print help message if an invalid string was used to specify
## the input or output file format
def help_wrong_file_format(
	fileFormat
    ):
    """
    Prints a message that the wrong file format was specified
    and then prints out a list of acceptable file formats

    Arguments
    ---------
    argv: fileFormat
        user specified input or output file format  
    """

    print("\n")
    print(fileFormat, "is not an accepted alignment format. Accepted alignment formats are:")
    print("'fasta', 'clustal', 'maf', 'mauve', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm'\n")
    print("For more help, use the -h/--help argument\n")
