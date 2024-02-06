import textwrap
import time
from .logger import logger
from .stats import TrimmingStats

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .files import FileFormat
    from .helpers import SeqType
    from .modes import TrimmingMode


def write_processing_aln() -> None:
    logger.info(
        textwrap.dedent(
            """\
        
        ------------------------
        | Processing Alignment |
        ------------------------
    """  # noqa
        )
    )


def write_user_args(
    in_file_name: str,
    in_file_format: "FileFormat",
    out_file_name: str,
    out_file_format: "FileFormat",
    seq_type: "SeqType",
    gaps: float,
    gap_chars: list,
    mode: "TrimmingMode",
    complement: bool,
    codon: bool,
    use_log: bool,
) -> None:
    if seq_type.value == "nt":
        seq_type_name = "Nucleotides"
    else:
        seq_type_name = "Protein"
    """
    Function to print user arguments to stdout
    """
    logger.info(
        textwrap.dedent(
            f"""\

    -------------
    | Arguments |
    -------------
    Input file: {in_file_name} (format: {in_file_format.value})
    Output file: {out_file_name} (format: {out_file_format.value})
    Sequence type: {seq_type_name}
    Gaps threshold: {gaps}
    Gap characters: {gap_chars}
    Trimming mode: {mode.value}
    Create complementary output: {complement}
    Process as codons: {codon}
    Create log file: {use_log}
    """  # noqa
        )
    )


def write_output_files_message(
    out_file_name: str, complement: bool, use_log: bool
) -> None:
    """
    Function to print out that the output files are being written
    """
    logger.info(
        textwrap.dedent(
            f"""\


        ------------------------
        | Writing output files |
        ------------------------
        Trimmed alignment: {out_file_name}
        Complement file: {out_file_name + '.complement' if complement else False}
        Log file: {out_file_name + '.log' if use_log else False}
    """
        )
    )


def write_output_stats(stats: "TrimmingStats", start_time: float) -> None:
    """
    Function to print out output statistics
    """
    logger.info(
        textwrap.dedent(
            f"""\

        ---------------------
        | Output Statistics |
        ---------------------
        Number of sites kept: {stats.output_length}
        Number of sites trimmed: {stats.trimmed_length}

        Percentage of alignment trimmed: {stats.trimmed_percentage}%

        Execution time: {round(time.time() - start_time, 3)}s
    """
        )
    )
