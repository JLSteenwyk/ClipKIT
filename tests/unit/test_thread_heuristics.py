from multiprocessing import cpu_count

from clipkit.clipkit import determine_effective_threads
from clipkit.modes import TrimmingMode


def test_effective_threads_respects_single_thread_request():
    effective = determine_effective_threads(
        requested_threads=1,
        mode=TrimmingMode.kpi,
        n_sequences=1500,
        alignment_length=30000,
    )
    assert effective == 1


def test_effective_threads_downshifts_kpi_small_workload():
    effective = determine_effective_threads(
        requested_threads=4,
        mode=TrimmingMode.kpi,
        n_sequences=1000,
        alignment_length=20000,
    )
    assert effective == 1


def test_effective_threads_limits_kpi_moderate_workload():
    effective = determine_effective_threads(
        requested_threads=8,
        mode=TrimmingMode.kpi,
        n_sequences=2000,
        alignment_length=60000,
    )
    assert effective == 2


def test_effective_threads_limits_kpi_large_workload():
    effective = determine_effective_threads(
        requested_threads=12,
        mode=TrimmingMode.kpic,
        n_sequences=3000,
        alignment_length=120000,
    )
    assert effective <= 4
    assert effective >= 1


def test_effective_threads_keeps_requested_for_non_kpi_modes():
    requested_threads = 4
    effective = determine_effective_threads(
        requested_threads=requested_threads,
        mode=TrimmingMode.gappy,
        n_sequences=1500,
        alignment_length=30000,
    )
    assert effective == min(requested_threads, max(1, cpu_count() or 1))
