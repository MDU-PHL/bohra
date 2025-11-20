import click
import pathlib
import os

from bohra.launcher.CheckDeps import _check_databases


@click.command()
def init_databases():
    """
    Check that dependencies are installed correctly.
    """
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    
    print("Will now check you installation and install anything that may be required. ")
    _check_databases(db_install=True)
    