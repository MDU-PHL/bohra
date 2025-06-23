"""Bohra 
.. moduleauthor:: Kristy Horan <kristyhoran15@gmail.com>

Bohra exinct tree kangaroo that lived on the nullarbor

Bohra is microbial genomics pipeline, designed predominantly for use in public health, but may also be useful in research settings. The pipeline takes as input a tab-delimited file with the isolate IDs followed by the path to READ1 and READ2, a reference for alignment and a unique identifier, where reads are illumina paired end reads (other platforms are not supported).

The pipline is based on nullarbor (https://github.com/tseemann/nullarbor) and is designed to be run in high performance computing environment.

Bohra is modular allowing the user to choose between calling SNPs and generating a phylogenetic tree, perform assemblies and detect AMR, perform typing etc; Or use the full pipeline to call SNPs, generate phylogenies, assemble, type and detect the pan-geneome. The output of Bohra is a html report that can be distributed, with downloadable tables and data.
"""


# import logging
# import argparse
import click

from bohra.SnpDetection import RunSnpDetection, SetupInputFiles, InitBohra, TestBohra
from bohra.version import version
from click.exceptions import UsageError
from click._compat import get_text_stderr

from bohra.launcher.Utils import _get_cmd_options
# from bohra.commands.pipelines.mash import mash
# from bohra.commands.pipelines.basic import basic
# from bohra.commands.pipelines.assemble import assemble
# from bohra.commands.pipelines.default import default
# from bohra.commands.pipelines.snps import snps
# from bohra.commands.pipelines.amr_typing import amr_typing
# from bohra.commands.pipelines.ska import ska
# from bohra.commands.pipelines.full import full
# from bohra.commands.pipelines.custom import custom
# from bohra.commands.pipelines.tb import tb
from bohra.commands.check import check
from bohra.commands.install import install 
from bohra.commands.generate_input import generate_input
from bohra.commands.test import bohratest


def _show_usage_error(self, file=None):
    if file is None:
        file = get_text_stderr()
    color = None
    if self.ctx is not None:
        color = self.ctx.color
        click.echo(self.ctx.get_help() + '\n', file=file, color=color)
    click.echo('Error: %s' % self.format_message(), file=file, color=color)

UsageError.show = _show_usage_error

@click.group()
def cli():
    pass

@cli.group()
def run():
    """
    Run the Bohra pipeline.
    """
    pass



def create_subcommand_with_options(name, options_dict):
    f"""Dynamically created a subcommand with options from a list."""

    @run.command(name=name, help = f"Help for the {name} pipeline.")
    def run_subcommand(**kwargs):
        # """A dynamically generated subcommand {name}."""
        click.echo(f"Running {name} subcommand with options:")
        
        for key, value in kwargs.items():
            click.echo(f"  {key}: {value}")

    # Add options dynamically
    # print(options_list)
    for opt in options_dict:
        click.option(
            f"--{opt['name']}",
            type=opt.get('type', None),
            default=opt.get('default', None),
            help=opt.get('help', ''),
            show_default= True,
            # is_flag= opt.get('is_flag', False),
        )(run_subcommand) # Apply the decorator to the function

    return run_subcommand

cmd_opts = _get_cmd_options()

for opt in cmd_opts:
    create_subcommand_with_options(name = opt, options_dict = cmd_opts[opt])

cli.add_command(check.check)
cli.add_command(install.install_deps)
cli.add_command(generate_input.generate_input)
cli.add_command(bohratest.test)



if __name__ == '__main__':
    cli()

