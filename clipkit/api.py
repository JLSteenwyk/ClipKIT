from typing import TYPE_CHECKING, Union
from tempfile import NamedTemporaryFile

from .clipkit import run
from .files import FileFormat
from .helpers import SeqType, write_msa
from .logger import logger
from .modes import TrimmingMode
from .plot_report import write_trim_plot_report

if TYPE_CHECKING:
    from .clipkit import TrimRun
    from .stats import TrimmingStats


ClipkitReturn = Union[tuple["TrimRun", "TrimmingStats"], tuple[str, "TrimmingStats"]]


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
    sequence_type: Union[SeqType, str, None] = None,
    codon: bool = False,
    ends_only=False,
    threads: int = 1,
    plot_trim_report_path: Union[str, None] = None,
) -> ClipkitReturn:
    """
    If input_file_path is given with no output_file_path -> Bio MSA (multiple sequence alignment object)
    If input_file_path is given and output_file_path is given -> write to output file
    If raw_alignment is given we write it to NamedTemporaryFile and then pass to execute
        * handles when output_file_path is given and also when not given
    """
    if threads < 1:
        raise ValueError("threads must be an integer >= 1")
    has_raw_alignment = raw_alignment is not None
    has_input_path = input_file_path is not None
    if has_raw_alignment == has_input_path:
        raise ValueError(
            "Provide exactly one of raw_alignment or input_file_path."
        )
    if has_raw_alignment and raw_alignment == "":
        raise ValueError("raw_alignment cannot be empty.")
    if has_input_path and input_file_path == "":
        raise ValueError("input_file_path cannot be empty.")

    original_logger_disabled = logger.disabled
    logger.disabled = True
    output_temp_file = None
    input_temp_file = None
    try:
        if raw_alignment:
            input_temp_file = NamedTemporaryFile()
            input_temp_file.write(bytes(raw_alignment, "utf-8"))
            input_temp_file.flush()

        if not output_file_path:
            output_temp_file = NamedTemporaryFile()

        if isinstance(sequence_type, str):
            try:
                sequence_type = SeqType(sequence_type.lower())
            except ValueError as exc:
                raise ValueError(
                    "sequence_type must be one of: 'aa', 'nt', SeqType.aa, SeqType.nt, or None."
                ) from exc

        # override options not currently available through programmatic interface
        complement = False
        use_log = False
        quiet = True
        auxiliary_file = None

        trim_run, stats = run(
            input_temp_file.name if input_temp_file else input_file_path,
            input_file_format,
            output_temp_file.name if output_temp_file else output_file_path,
            output_file_format,
            auxiliary_file,
            sequence_type,
            gaps,
            gap_characters,
            complement,
            codon,
            TrimmingMode(mode),
            use_log,
            quiet,
            ends_only,
            threads,
        )

        if plot_trim_report_path:
            write_trim_plot_report(
                plot_trim_report_path,
                trim_run.msa,
                mode=TrimmingMode(mode).value,
                gaps=trim_run.gaps,
                sequence_type=trim_run.sequence_type.value,
            )

        if not output_file_path:
            return trim_run, stats

        base_metadata = trim_run.alignment.annotations.get("ecomp_metadata")
        write_msa(
            trim_run.msa,
            output_file_path,
            trim_run.output_file_format,
            base_metadata=base_metadata,
        )
        return output_file_path, stats
    finally:
        logger.disabled = original_logger_disabled
        if input_temp_file is not None:
            input_temp_file.close()
        if output_temp_file is not None:
            output_temp_file.close()
