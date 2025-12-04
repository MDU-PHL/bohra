import click
import pathlib
import os

from bohra.launcher.SetupInput import find_data

@click.command()
@click.option('--reads',
              default="",
              help="Path to search for reads files, e.g. *.f*q.gz",)
@click.option('--contigs',
              default="",
              help="Path to search for assembly files, e.g. *.f*a.gz")
@click.option('--isolate_ids',
              default='',
              help="Path to a file containing at least one column 'Isolate' with isolate names. Optionally add 'species' and other columns you wish to use for further annotation of trees.")
@click.option('--outname',
              default='bohra_input.tsv',
              help="Name of the file to write the generated input table to.")


def generate_input(reads, contigs, isolate_ids, outname):
    """
    Generare input files for the Bohra pipeline.
    """
    
    if not reads and not contigs:
        click.echo("Please specify  --reads and/or --contigs.")
        return

    
    # click.echo(f"Generating input files for {input_type} in {path}.")

    find_data(
        isolate_ids=isolate_ids,
        reads=reads,
        contigs=contigs,
        outname=outname
    )

  