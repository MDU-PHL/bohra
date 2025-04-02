import click
import pathlib
import os


@click.command()
@click.option('--reads',
              is_flag=True,
              help="Set if your input type is reads.")
@click.option('--contigs',
              is_flag=True,
              help="Set if your input type is contigs.")
@click.option('--path',
              default=pathlib.Path.cwd().absolute(),
              help='The directory where your input files are located, default is current directory',
              type=click.Path(exists=True))
@click.option('--isolate_ids',
              default='',
              help='File containing isolate IDs one per line - OPTIONAL.')
def generate_input():
    """
    Generare input files for the Bohra pipeline.
    """
    print("Generating input files...")
  