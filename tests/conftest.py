# global fixtures can go here
import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
    config.addinivalue_line("markers", "slow: mark as slow test")
