import click
import pathlib
import os


@click.command()
@click.option('--reads', '-r',
              help='Path to reads file, which is a tab-delimited with 3 columns <isolatename>  <path_to_read1> <path_to_read2>.',
              default='')
@click.option('--kraken_db', '-k',
              default=os.getenv("KRAKEN2_DEFAULT_DB", ''),
              metavar='KRAKEN2_DEFAULT_DB',
              show_default=True,
              help="Path to DB for use with kraken2, if no DB present speciation will not be performed.")
@click.option('--cpus',
              help='Number of max CPU cores to run, will define how many rules are run at a time, if 0 then the avail cpus will be determined at time of launch', 
              default=0)
@click.option('--minmap', '-mp',
              help='Snippy - minimum read mapping quality to consider.', 
              default='60')
@click.option('--basequal', '-bq',
              help='Snippy - Minimum base quality to consider.', 
              default='13')
@click.option('--minqual', '-mq',
              help='Snippy - minumum QUALITY in VCF column 6', 
              default='100')
@click.option('--minfrac', '-mf',
              help='Snippy - minumum proportion for variant evidence ', 
              default='0')
@click.option('--mincov', '-mc',
              help='Snippy - minimum site depth to for calling alleles.', 
              default='10')
@click.option('--cluster',
              is_flag=True, 
              help='Set if you want to do heirarchical clustering.')
@click.option('--cluster_method', '-cm',
              help='The clustering method to use, default is \'single-linkage\'', 
              type=click.Choice(['single-linkage', 'average', 'complete', 'centroid', 'median', 'ward', 'weighted']),
              default='single-linkage')
@click.option('--cluster_threshold', '-ct',
              help='Comma separated list of thresholds to use for clustering, default is \'10\'',
              type=str, 
              default='5,12,20')
@click.option('--tree_builder', '-tb',
              default='fasttree',
              help='The tree builder to use, default is \'fasttree\'',
              type=click.Choice(['fasttree', 'iqtree']))
@click.option('--workdir', '-w',
              default=pathlib.Path.cwd().absolute(), 
              help='The directory where Bohra will be run, default is current directory', 
              type=click.Path(exists=True))
@click.option('--no_phylo',
              is_flag=True, 
              help='Set if you do NOT want to generate a phylogentic tree.')
@click.option('--conda_path',       
              default=pathlib.Path(os.getenv('CONDA_PREFIX', '')), 
              help='The path to where your pre-installed conda envs are stored, defaults to installing conda envs in your work directory. This can be provided in your profiles settings as well - it assumes you have pre-configured all of your conda environments for each process run by bohra, this is an advanced setting. Please take care.')
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
def tb():
    """
    Run a Mtb pipeline using Bohra. This pipeline will run a series of steps to generate a phylogenetic tree and report for the provided reads.
    """
    print("Running bohra tb...")