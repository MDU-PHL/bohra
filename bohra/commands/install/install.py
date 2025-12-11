import click
import pathlib
import os

from bohra.launcher.CheckDeps import check_dependencies, _check_databases


@click.command()
@click.option('--force_reinstall',
              is_flag=True,
            #   default=False,
              help="Force reinstallation of all dependencies, even if they are already installed.")

@click.option('--install_dbs',
              is_flag=True,
            #   default=False,
              help="Install required databases.")
@click.option('--setup_databases',
              is_flag=True,
            #   default=False,
              help="Download and/or setup required databases.")
@click.option('--tool',
              # is_flag=True,
              type=click.Choice(['all', 'torstyverse', 'seqquality', 'relationships','snippy', 'ectyper','mob_suite','panaroo','kleborate','stype','tamr','sonneitype','classify-pangenome'], case_sensitive=False),
              default="all",
              help="Install only a specific set of tools from a single environment. Should really only be used for development and/or testing purposes.")

def install_deps(force_reinstall:bool=False, install_dbs:bool=False, setup_databases:bool=False, tool:str="all"):
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    print("Will now check and install anything that may be required. ")
    # print(str(force_reinstall).lower())
    check_dependencies(force_reinstall=str(force_reinstall).lower(), tool =tool )
    if install_dbs or setup_databases:
        print(f"Will now check your databases{' and install anything that may be required' if setup_databases else ''}. ")
        _check_databases(db_install=setup_databases)
