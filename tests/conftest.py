# global fixtures can go here
import pytest
import os
from pathlib import Path


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
    config.addinivalue_line("markers", "slow: mark as slow test")


@pytest.fixture(autouse=True)
def ensure_output_directory():
    """Ensure output directory exists for tests"""
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)
    yield
    # Optionally clean up after tests (commented out to preserve test outputs)
    # if output_dir.exists():
    #     import shutil
    #     shutil.rmtree(output_dir)
