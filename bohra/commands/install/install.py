import click
import pathlib
import os

from bohra.launcher.CheckDeps import check_dependencies


@click.command()

def install_deps():
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    print("Will now check and install anything that may be required. ")
    check_dependencies( )
