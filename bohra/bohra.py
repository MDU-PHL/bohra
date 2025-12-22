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
import sys

from bohra.version import version
from click.exceptions import UsageError
from click._compat import get_text_stderr

from bohra.launcher.Utils import _get_run_cmd_options, _get_dep_cmd_options
from bohra.commands.init_databases import init_databases
# from bohra.commands.install import install 
from bohra.commands.generate_input import generate_input
from bohra.commands.test import bohratest
# from bohra.commands.check_setup import check_setup
from bohra.launcher.BohraRun import run_bohra
from bohra.launcher.Deps import dependencies


def _show_usage_error(self, file=None):
    if file is None:
        file = get_text_stderr()
    color = None
    if self.ctx is not None:
        color = self.ctx.color
        click.echo(self.ctx.get_help() + '\n', file=file, color=color)
    click.echo('ERROR: %s' % self.format_message(), file=file, color=color)

UsageError.show = _show_usage_error

@click.group(no_args_is_help=True,invoke_without_command=True)
@click.version_option(version, message='%(prog)s version %(version)s')
# @click.pass_context
def cli():
    pass

@cli.group(no_args_is_help=True,invoke_without_command=True)
def run():
    """
    Run the Bohra pipeline.
    """
    pass


@cli.group(no_args_is_help=True,invoke_without_command=True)
def deps():
    """
    Manage bohra dependencies.
    """
    pass


def generic_options(_func: click.Command, options: dict) -> click.Command:
    """Add generic options to a Click command."""
    for opt in options:
        if 'short_name' in opt:
            click.option(
                f"--{opt['name']}",
                f"{opt['short_name']}",
                type=opt.get('type', None),
                default=opt.get('default', None),
                help=opt.get('help', ''),
                show_default= True,
                is_flag= opt.get('is_flag', False),
            )(_func) # Apply the decorator to the function
        else:
            click.option(
                f"--{opt['name']}",
                type=opt.get('type', None),
                default=opt.get('default', None),
                help=opt.get('help', ''),
                show_default= True,
                is_flag= opt.get('is_flag', False),
            )(_func)
    return _func
    

def create_deps_subcommand_with_options(name, options_dict):
    f"""Dynamically created a subcommand with options from a list."""

    @deps.command(name=name, help = f"Help for {name}ing dependencies.")
    def deps_subcommand(**kwargs):
        try:
            if dependencies(_action=name, kwargs=kwargs) == 0:
                click.echo(f"Successfully completed {name}ing dependencies.")
            else:
                raise UsageError(f"An error occurred while {name}ing dependencies.")
        except Exception as e:
            raise UsageError(f"An error occurred while installing the {name} dependencies: {e}")


    deps_subcommand = generic_options(deps_subcommand, options_dict)
    return deps_subcommand


def create_run_subcommand_with_options(name, options_dict):
    f"""Dynamically created a subcommand with options from a list."""

    @run.command(name=name, help = f"Help for the {name} pipeline.")
    def run_subcommand(**kwargs):
        # try:
        if len(sys.argv) > 3:
            run_bohra(pipeline=name, kwargs=kwargs)
        else:
            raise UsageError(f"You did not supply any input. Please supply some input to run the {name} pipeline.")
        # except Exception as e:
        #     raise UsageError(f"An error occurred while running the {name} pipeline: {e}")

    # for opt in options_dict:
        
    #     if 'short_name' in opt:
    #         click.option(
    #             f"--{opt['name']}",
    #             f"{opt['short_name']}",
    #             type=opt.get('type', None),
    #             default=opt.get('default', None),
    #             help=opt.get('help', ''),
    #             show_default= True,
    #             is_flag= opt.get('is_flag', False),
    #         )(run_subcommand) # Apply the decorator to the function
    #     else:
    #         click.option(
    #             f"--{opt['name']}",
    #             type=opt.get('type', None),
    #             default=opt.get('default', None),
    #             help=opt.get('help', ''),
    #             show_default= True,
    #             is_flag= opt.get('is_flag', False),
    #         )(run_subcommand)
    run_subcommand = generic_options(run_subcommand, options_dict)
    return run_subcommand

run_cmd_opts = _get_run_cmd_options()
dep_cmd_opt = _get_dep_cmd_options()
for opt in run_cmd_opts:
    create_run_subcommand_with_options(name = opt, options_dict = run_cmd_opts[opt])

for opt in dep_cmd_opt:
    create_deps_subcommand_with_options(name = opt, options_dict = dep_cmd_opt[opt])

cli.add_command(init_databases.init_databases)
# cli.add_command(install.install_deps)
# cli.add_command(check_setup.check_deps)
cli.add_command(generate_input.generate_input)
cli.add_command(bohratest.test)



if __name__ == '__main__':
    cli()

