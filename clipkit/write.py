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
    inFile, inFileFormat, outFile, outFileFormat, char, gaps, mode, complement, use_log
):
    if char.value == "nt":
        seq_type = "Nucleotides"
    else:
        seq_type = "Protein"
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
    Sequence type: {seq_type}
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


def write_output_stats(stats: dict, start_time):
    """
    Function to print out output statistics
    """
    print(
        textwrap.dedent(
            f"""\

        ---------------------
        | Output Statistics |
        ---------------------
        Number of sites kept: {stats["output_length"]}
        Number of sites trimmed: {stats["trimmed_length"]}

        Percentage of alignment trimmed: {stats["trimmed_percentage"]}%

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
