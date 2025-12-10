import click
import pathlib
import os
from bohra.launcher.TestBohra import run_tests

@click.command()
@click.option('--cpus',
              default=1,
              help="Number of CPUs to use for testing Bohra installation.") 
def test(cpus:int=1):
    """
    Check that bohra is installed correctly and runs as expected.
    """
    
    run_tests(cpus=cpus)