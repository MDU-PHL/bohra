import click


from bohra.launcher.InstallDeps import install_dependencies

@click.command()
@click.option('--prefix','-p',
              default='bohra',
              help='The prefix for your environments, this will be used to create the conda environments for each process. The default is \'bohra\'')
@click.option('--install-deps/--no-install-deps', 
              is_flag=True, 
              default=True,
              help = "Creates fixed conda environments for each process. This is highly recommended. If you do not want to install dependencies, you can use --no-install-deps.")
@click.option('--databases/--no-databases',
              is_flag=True,
              default=True,
              help = "Checks for databases required for the pipeline and that environment variables are set properly. This is highly recommended. If you do not want to check for databases, you can use --no-databases.")
def install_deps(prefix, install_deps, databases):
    """
    Install dependencies for Bohra - Highly recommended to run this before running the pipeline.
    """
    print("Will now install dependencies for Bohra.")
    print("This will take some time.")
    print("Please be patient.")
    if install_deps:
        print("Installing dependencies...")
        install_dependencies(prefix=prefix)
    else:
        print("Skipping installation of dependencies.")
    if databases:
        print("Checking for databases...")
    else:
        print("Skipping database checks.")  
    # Add the code to check dependencies here
    # This is a placeholder for the actual implementation
    