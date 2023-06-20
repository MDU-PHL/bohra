import pytest

def pytest_configure(config):
    config.addinivalue_line(
        "markers", "nf: marks tests as nextflow (deselect with '-m \"not nextflow\"')"
    )