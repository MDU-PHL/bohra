import click
import pathlib
import os

@click.command()
# common to all pipelines
@click.option('--input_file', '-i',
              help='Path to reads file, which is a tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>.',
              default='')
@click.option('--sylph_db', '-k',
              default="",
              show_default=True,
              help="Path to DB for use with sylph")
@click.option('--kraken_db', '-k',
              default=os.getenv("KRAKEN2_DEFAULT_DB", ''),
              metavar='KRAKEN2_DEFAULT_DB',
              show_default=True,
              help="Path to DB for use with kraken2")
@click.option('--speciation/--no-speciation',
              is_flag=True, 
              help='Speciation will be performed by deafult - use --no-speciation if you do not need species detected.')

# resource options common to all pipelines
@click.option('--cpus',
              help='Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch', 
              default=0)
@click.option('--workdir', '-w',
              default=pathlib.Path.cwd().absolute(), 
              help='The directory where Bohra will be run, default is current directory', 
              type=click.Path(exists=True))
@click.option('--conda_path',       
              default=pathlib.Path(os.getenv('CONDA_PREFIX', '')), 
              help='The path to where your pre-installed conda bohra-envs are stored. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care if you are changing it.')
@click.option('--keep',
              default='N', 
              type=click.Choice(['Y', 'N']), 
              help='If you are rerunning bohra over an exisiting directory set --keep to \'Y\' to archive report files - otherwise previous report files will be removed.')
@click.option('--proceed',
              is_flag=True, 
              help='If you would like to proceed straigt to the pipeline.')
@click.option('--force', '-f',
              is_flag=True, 
              help='Add if you would like to force a complete restart of the pipeline. All previous logs will be lost.')
@click.option('--no-conda',
              is_flag=True, 
              help='Set if you DO NOT WANT to use separate conda environments for each nextflow process.')
@click.option('--nfconfig','-nfcfg',
              default = f"", 
              help='An additional config file, required if running on a non-local machine, ie slurm, cloud. For help see documentation at https://github.com/MDU-PHL/bohra or https://www.nextflow.io/docs/latest/executor.html',) # don't need this
@click.option('--profile',
              default=f"", 
              help='The resource profile to use. Defaults to local, if using an alternative config file, this value should represent the name of a profile provided')
def basic(input_file: str,
          sylph_db: str,
          kraken_db: str,
          speciation: bool,
          conda_path: str,
          keep: str,
          proceed: bool,
          force: bool,
          no_conda: bool,
          nfconfig: str,
          profile: str):
    """
    Run bohra for basic sequence assessment
    """
    print("Running basic...")