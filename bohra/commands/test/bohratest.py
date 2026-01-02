import click
import pathlib
import os
from bohra.launcher.TestBohra import run_tests

@click.command()
@click.option('--cpus',
              default=1,
              help="Number of CPUs to use for testing Bohra installation.") 
@click.option('--shovill_ram',
              default=16,
              help="Amount of RAM to allocate to shovill assembler.") 
@click.option('--wdir',
              default=str(pathlib.Path.cwd()),
              help="Working directory for the test run. Default is the current working directory.")
def test(cpus:int=1, shovill_ram:int=16, wdir:str=str(pathlib.Path.cwd())):
    """
    Check that bohra is installed correctly and runs as expected.
    """
    
    run_tests(cpus=cpus, shovill_ram=shovill_ram, wdir=wdir)