import click
import pathlib
import os

from bohra.launcher.CheckDeps import install_dependencies, check_databases


@click.command()
@click.option('--prefix','-p',
              default=f"{pathlib.Path(os.environ['CONDA_PREFIX']).name}" if os.environ['CONDA_PREFIX'] != "" else 'bohra',
              help='The prefix for your environments, this will be used to create the conda environments for each process.',
              show_default=True)
@click.option('--install-deps/--no-install-deps', 
              is_flag=True, 
              default=True,
              help = "Creates fixed conda environments for each process. This is highly recommended. If you do not want to install dependencies, you can use --no-install-deps.")
@click.option('--databases/--no-databases',
              is_flag=True,
              default=True,
              help = "Checks for databases required for the pipeline and that environment variables are set properly. This is highly recommended. If you do not want to check for databases, you can use --no-databases.")

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
    if install_deps:
        print("Will now check you installation and install anything that may be required. ")
        install_dependencies(prefix=prefix, check_only = False)
    else:
        print("Will just check your installation - no dependencies will be installed.")
        install_dependencies(prefix=prefix, check_only = True)
    if databases:
        print("Checking for databases...")
        check_databases(prefix=prefix)
    else:
        print("Skipping database checks.") 