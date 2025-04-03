import click
import pathlib
import os
from bohra.launcher.TestBohra import run_tests

@click.command()
def test():
    """
    Check that bohra is installed correctly and runs as expected.
    """
    
    run_tests()