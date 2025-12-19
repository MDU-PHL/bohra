import click
import pathlib
import os

from bohra.launcher.Deps import _check_databases


@click.command()
@click.option('--setup_databases',
              is_flag=True,
            #   default=False,
              help="Download and/or setup required databases.")
def init_databases(setup_databases:bool=False):
    """
    Check that dependencies are installed correctly.
    """
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    db_install = True if setup_databases else False
    print(f"Will now check your databases{' and install anything that may be required' if db_install else ''}. ")
    _check_databases(db_install=db_install)
    