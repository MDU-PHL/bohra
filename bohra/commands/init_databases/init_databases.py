import click
import pathlib
import os

from bohra.launcher.Deps import _check_databases


@click.command()
@click.option('--setup_databases',
              is_flag=True,
            #   default=False,
              help="Download and/or setup required databases.")
@click.option('--database_path',
              default="",
              help="If you select to setup databases, specify the path to download them to.")
def init_databases(setup_databases:bool=False, database_path:str=""):
    
    db_install = True if setup_databases else False
    print(f"Will now check your databases{' and install anything that may be required' if db_install else ''}. ")
    _check_databases(db_install=db_install, database_path=database_path)
    