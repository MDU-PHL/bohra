import click
import pathlib
import os

from bohra.launcher.SetupInput import find_data

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

def generate_input(reads, contigs, path, isolate_ids):
    """
    Generare input files for the Bohra pipeline.
    """
    if reads and contigs:
        click.echo("Please specify either --reads or --contigs, not both.")
        return

    if not reads and not contigs:
        click.echo("Please specify either --reads or --contigs.")
        return

    if reads:
        input_type = 'reads'
    else:
        input_type = 'contigs'

    # Check if the path exists
    if not os.path.exists(path):
        click.echo(f"Path {path} does not exist.")
        return

    # Check if the isolate_ids file exists
    if isolate_ids and not os.path.exists(isolate_ids):
        click.echo(f"File {isolate_ids} does not exist.")
        return

    # click.echo(f"Generating input files for {input_type} in {path}.")

    find_data(
        input_type=input_type,
        isolate_ids=isolate_ids,
        path=path
    )

  