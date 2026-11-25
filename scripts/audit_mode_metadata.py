"""Emit exact trimming metadata for the selected ClipKIT source tree."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from clipkit import __file__ as clipkit_path
from clipkit.clipkit import run
from clipkit.helpers import SeqType
from clipkit.modes import AmbiguityHandling, TrimmingMode
from run_benchmark_smoke import (
    CASES,
    FILES,
    _array_digest,
    _hash_bytes,
    _write_generated_fixtures,
)

AUDIT_CASES = (
    "gappy_small",
    "smart_gap_medium",
    "block_gappy_small",
    "gappyout_small",
    "kpi_small",
    "kpic_small",
    "kpi_gappy_small",
    "kpic_gappy_small",
    "kpi_smart_gap_small",
    "kpic_smart_gap_small",
    "entropy_small",
    "composition_bias_small",
    "heterotachy_medium",
    "cst_tiny",
    "c3_medium_nt",
    "cli_entropy_missing_medium_nt",
    "cli_entropy_fractional_medium_nt",
    "cli_entropy_literal_medium_nt",
    "codon_medium_nt",
    "ends_only_medium_nt",
    "stop_terminal_end_to_end",
    "stop_internal_end_to_end",
    "stop_all_end_to_end",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path, help="Source tree added to PYTHONPATH")
    source = parser.parse_args().source_root.resolve()
    if Path(clipkit_path).resolve().parent.parent != source:
        parser.error("PYTHONPATH must select the requested ClipKIT source tree")

    results = {}
    with TemporaryDirectory(prefix="clipkit-metadata-audit-") as directory:
        generated = _write_generated_fixtures(Path(directory))
        for name in AUDIT_CASES:
            case = CASES[name]
            dataset = FILES[case.dataset]
            path = source / dataset if dataset else generated[case.dataset]
            auxiliary = str(source / case.auxiliary_file) if case.auxiliary_file else None
            result, stats = run(
                input_file=str(path),
                input_file_format=None,
                output_file=None,
                output_file_format=None,
                auxiliary_file=auxiliary,
                sequence_type=SeqType(case.sequence_type) if case.sequence_type else None,
                gaps=case.gaps,
                gap_characters=None,
                complement=False,
                codon=case.codon,
                mode=TrimmingMode[case.mode],
                use_log=False,
                quiet=True,
                ends_only=case.ends_only,
                threads=case.threads,
                remove_stop_codons=case.stop_codons,
                ambiguity_handling=AmbiguityHandling(case.ambiguity_handling or "missing"),
            )
            msa = result.msa
            records = [
                (record.id, record.name, record.description, str(record.seq))
                for record in result.trimmed
            ]
            classifications = "\n".join(
                value.value for value in msa.site_classification_types
            )
            results[name] = {
                "mode": case.mode,
                "dataset": case.dataset,
                "threshold": float(result.gaps),
                "stats": stats.summary,
                "stop_codon_masking": result.stop_codon_masking.summary,
                "output_sha256": _hash_bytes(
                    json.dumps(records, ensure_ascii=False).encode("utf-8")
                ),
                "keep_sha256": _array_digest(msa._site_positions_to_keep, "int64"),
                "trim_sha256": _array_digest(msa._site_positions_to_trim, "int64"),
                "classification_sha256": _hash_bytes(classifications.encode("utf-8")),
                "gappyness_sha256": _array_digest(msa.site_gappyness),
            }

    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
