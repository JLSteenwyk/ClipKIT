import runpy
from pathlib import Path

import pytest

from clipkit.version import __version__


@pytest.mark.integration
def test_repository_api_example_remains_runnable(capsys):
    script_path = Path(__file__).parents[2] / "test-api.py"

    runpy.run_path(str(script_path), run_name="clipkit_api_example")

    assert capsys.readouterr().out.strip() == __version__
