import click
import pathlib
import os

from bohra.launcher.CheckDeps import check_dependencies, _check_databases


@click.command()

@click.option('--install_dbs',
              is_flag=True,
            #   default=False,
              help="Force reinstallation of all dependencies, even if they are already installed.")
@click.option('--setup_databases',
              is_flag=True,
            #   default=False,
              help="Download and/or setup required databases.")

def check_deps(install_dbs:bool=False, setup_databases:bool=False):
    """
    Check dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now check dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    print("Will now check and install anything that may be required. ")
    # print(str(force_reinstall).lower())
    check_dependencies(force_reinstall='false' )
    _check_databases(db_install=install_dbs or setup_databases)
