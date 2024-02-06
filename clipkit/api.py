from typing import TextIO, Union
from tempfile import NamedTemporaryFile

from .clipkit import run
from .files import FileFormat
from .helpers import SeqType, write_msa
from .logger import logger
from .modes import TrimmingMode


def clipkit(
    *,
    raw_alignment: Union[str, None] = None,
    input_file_path: Union[str, None] = None,
    output_file_path: Union[str, None] = None,
    mode: TrimmingMode = TrimmingMode.smart_gap,
    gaps: Union[float, None] = None,
    gap_characters=None,
    input_file_format=FileFormat.fasta,
    output_file_format=FileFormat.fasta,
    sequence_type=SeqType.aa
) -> TextIO:
    """
    If input_file_path is given with no output_file_path -> Bio MSA (multiple sequence alignment object)
    If input_file_path is given and output_file_path is given -> write to output file
    If raw_alignment is given we write it to NamedTemporaryFile and then pass to execute
        * handles when output_file_path is given and also when not given
    """
    logger.disabled = True
    output_temp_file = None
    input_temp_file = None
    if raw_alignment:
        input_temp_file = NamedTemporaryFile()
        input_temp_file.write(bytes(raw_alignment, "utf-8"))
        input_temp_file.flush()

    if not output_file_path:
        output_temp_file = NamedTemporaryFile()

    # override some options not currently available through programmatic interface
    complement = False
    codon = False
    use_log = False
    quiet = True

    trim_run, stats = run(
        input_temp_file.name if input_temp_file else input_file_path,
        input_file_format,
        output_temp_file.name if output_temp_file else output_file_path,
        output_file_format,
        sequence_type,
        gaps,
        gap_characters,
        complement,
        codon,
        TrimmingMode(mode),
        use_log,
        quiet,
    )

    if not output_file_path:
        return trim_run, stats
    else:
        write_msa(trim_run.msa, output_file_path, trim_run.output_file_format)
        return output_file_path, stats
