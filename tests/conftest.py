import pytest

def pytest_configure(config):
    config.addinivalue_line(
        "markers", "nf_typing: marks tests as nextflow typing modules (deselect with '-m \"not nextflow\"')"
    )
    config.addinivalue_line(
        "markers", "nf_snps: marks tests as nextflow snps modules m(deselect with '-m \"not nextflow\"')"
    )
    config.addinivalue_line(
        "markers", "setup_bohra: marks tests as setup tests m(deselect with '-m \"not nextflow\"')"
    )