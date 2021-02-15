import textwrap
import time


def write_processing_aln():
    """
    Function to print out processing alignment
    """
    print(
        textwrap.dedent(
            f"""\
        
        ------------------------
        | Processing Alignment |
        ------------------------
    """
        )
    )


def write_user_args(
    inFile, inFileFormat, outFile, outFileFormat, gaps, mode, complement, use_log
):
    """
    Function to print user arguments to stdout
    """
    print(
        textwrap.dedent(
            f"""\
    -------------
    | Arguments |
    -------------
    Input file: {inFile} (format: {inFileFormat.value})
    Output file: {outFile} (format: {outFileFormat.value})
    Gaps threshold: {gaps}
    Trimming mode: {mode.value}
    Create complementary output: {complement}
    Create log file: {use_log}
    """
        )
    )


def write_output_files_message(outFile, complement, use_log):
    """
    Function to print out that the output files are being written
    """
    print(
        textwrap.dedent(
            f"""\


        ------------------------
        | Writing output files |
        ------------------------
        trimmed alignment: {outFile}
        complement file: {outFile + '.complement' if complement else False}
        log file: {outFile + '.log' if use_log else False}
    """
        )
    )


def write_output_stats(alignment, keepD, trimD, start_time):
    """
    Function to print out output statistics
    """
    alignment_length = alignment.get_alignment_length()
    output_len = len(next(iter(keepD.values())))
    trimmed_len = len(next(iter(trimD.values())))
    print(
        textwrap.dedent(
            f"""\

        ---------------------
        | Output Statistics |
        ---------------------
        Number of sites kept: {output_len}
        Number of sites trimmed: {trimmed_len}

        Percentage of alignment trimmed: {round((trimmed_len / alignment_length) * 100, 3)}%

        Execution time: {round(time.time() - start_time, 3)}s
    """
        )
    )

def write_determining_smart_gap_threshold():
    """
    Function to inform using that the smart-gap threshold is being determined
    """
    print(
        textwrap.dedent(
            f"""\

        Determining smart-gap threshold...
    """
        )
    )
