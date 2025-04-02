import click
import pathlib
import os


@click.command()
@click.option('--prefix','-p',
              default='bohra',
              help='The prefix for your environments, this will be used to create the conda environments for each process. The default is \'bohra\'')
def install_deps():
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    # Add the code to check dependencies here
    # This is a placeholder for the actual implementation
    