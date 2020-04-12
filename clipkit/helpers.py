from Bio import AlignIO

####################################################################
### Supporting Functions 	                                     ###
### This block of code contains all supporting functions that    ###
### read the input files, calculate gappyness, etc  			 ###
####################################################################
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
   

## Function to get the sequence at a given column. Function
## will also determine features of the position including 
## length and number of gaps.
## Support function for keep_trim_and_log()
def get_sequence_at_position_and_report_features(
    alignment, i
    ):
    """
    Count the occurence of each character at a given position
    in an alignment. This information is used to determine
    if the aligment is parsimony informative or not. When 
    count characters, gaps ('-') are excluded

    Arguments
    ---------
    argv: seqAtPosition
        string that contains the sequence at a given column 
    """

    # save the sequence at the position to a string 
    seqAtPosition  = ''
    seqAtPosition += alignment[:, i]
    seqAtPosition  = seqAtPosition.upper()

    # determine the length and number of gaps in an alignment position
    lengthOfSeq = len(seqAtPosition)
    numOfGaps   = seqAtPosition.count('-')
    gappyness   = numOfGaps/lengthOfSeq

    return seqAtPosition, gappyness


## Function to count the number of occurences in each character
## Support function for keep_trim_and_log()
def count_characters_at_position(
    seqAtPosition
    ):
    """
    Count the occurence of each character at a given position
    in an alignment. This information is used to determine
    if the aligment is parsimony informative or not. When 
    count characters, gaps ('-') are excluded

    Arguments
    ---------
    argv: seqAtPosition
        string that contains the sequence at a given column 
    """

    numOccurences = {}
    for char in set(seqAtPosition.replace('-','')): 
        numOccurences[char] = seqAtPosition.count(char)
    return numOccurences
    

## Function to determine which positions of an alignment should be 
## kept or trimmed 
def keep_trim_and_log(
    alignment, gaps, mode
    ):
    """
    Determines positions to keep or trim and saves these positions
    to dictionaries named keepD and trimD. For both dictionaries,
    keys are the names of the sequence entries and the values are
    sequences that will be kept or trimmed

    Arguments
    ---------
    argv: alignment
        biopython multiple sequence alignment object
    argv: gaps
        gaps threshold to determine if a position is too gappy or not
    """

    # initialize dictionaries that will eventually be populated with
    # alignment positions to keep or trim (keys) and the sequence at
    # that position (values). Also, initialize a list of log information
    # that will be kept in an array format
    keepD = {}
    for entry in alignment:
        keepD[entry.id] = []
    trimD = {}
    for entry in alignment:
        trimD[entry.id] = []
    logArr = []

    # loop through alignment
    for i in range(0, alignment.get_alignment_length(), int(1)):
        
        # save the sequence at the position to a string and calculate the gappyness of the site
        seqAtPosition, gappyness = get_sequence_at_position_and_report_features(alignment, i)

        ## determine if the site is parsimony informative 
        ## "A site is parsimony-informative if it contains at least two types of nucleotides 
        ## (or amino acids), and at least two of them occur with a minimum frequency of two."
        ## https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

        # Create a dictionary that tracks the number of occurences of each character 
        # excluding gaps or '-'
        numOccurences = count_characters_at_position(seqAtPosition)

        # if the number of values that are greater than two 
        # in the numOccurences dictionary is greater than two, 
        # the site is parsimony informative
        d = dict((k, v) for k, v in numOccurences.items() if v >= 2)
        if len(d) >= 2:
            parsimony_informative = True
        else:
            parsimony_informative = False

        if mode == 'kpi-gappy':
            # if gappyness is lower than gaps threshold and the site is
            # parismony informative, save to keepD. Otherwise, save to trimD.
            # All the while, save information to logD
            if gappyness <= gaps and parsimony_informative:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                temp = []
                temp.append(str(i+1))
                temp.append('keep')
                temp.append("PI")
                temp.append(str(gappyness))
                logArr.append(temp)
            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                temp = []
                temp.append(str(i+1))
                temp.append('trim')
                temp.append("nPI")
                temp.append(str(gappyness))
                logArr.append(temp)
        elif mode == 'gappy':
            # if gappyness is lower than gaps threshold and the site is
            # parismony informative, save to keepD. Otherwise, save to trimD.
            # All the while, save information to logD
            if gappyness <= gaps:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                if parsimony_informative:
                    temp = []
                    temp.append(str(i+1))
                    temp.append('keep')
                    temp.append("PI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
                else: 
                    temp = []
                    temp.append(str(i+1))
                    temp.append('keep')
                    temp.append("nPI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                if parsimony_informative:
                    temp = []
                    temp.append(str(i+1))
                    temp.append('trim')
                    temp.append("PI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
                else: 
                    temp = []
                    temp.append(str(i+1))
                    temp.append('trim')
                    temp.append("nPI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
        
        elif mode == 'kpi':
            # if site is parismony informative, save to keepD. 
            # Otherwise, save to trimD.
            # All the while, save information to logD
            if parsimony_informative:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
 
                temp = []
                temp.append(str(i+1))
                temp.append('keep')
                temp.append("PI")
                temp.append(str(gappyness))
                logArr.append(temp)

            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])

                temp = []
                temp.append(str(i+1))
                temp.append('trim')
                temp.append("nPI")
                temp.append(str(gappyness))
                logArr.append(temp)
            
    # join elements in value lists in keepD and trimD 
    for k, v in keepD.items():
        keepD[k] = ''.join(v)
    for k, v in trimD.items():
        trimD[k] = ''.join(v)

    return(keepD, trimD, logArr)