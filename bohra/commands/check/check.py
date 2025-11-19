import click
import pathlib
import os

from bohra.launcher.CheckDeps import check_dependencies


@click.option('--envs',
              default=f"{pathlib.Path(__file__).parent.parent}/environments",
              help='The prefix for your environments, this will be used to create the conda environments for each process.',
              )

@click.command()
def check(prefix, install_deps, databases):
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
    check_dependencies()
    