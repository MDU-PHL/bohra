#!/usr/bin/env bash

#set -e -u -o pipefail

# check parameters
if [ "$#" -ne 3 ]; then
	EXE=$(basename $0)
	echo "USAGE: $EXE <envsdirr=path> <action=install|check> <force=true|false>"
	exit 255
fi

ENVSDIR=$1
ACTION=$2
FORCE=$3

BOHRA_CONDA_ENVS="$CONDA_PREFIX/conda_envs"
#BOHRA_CONDA_ENVS="$HOME/.conda/envs"
INSTALLER="conda"

# resets to base env
eval "$(conda shell.bash hook)"
$INSTALLER -V

# check that mamba is installed
if [ -x "$(which mamba)" ] ; then
    INSTALLER="mamba"
fi
echo "Using installer: $INSTALLER"

# print a bold string
function print_bold {
    echo -e "\033[1m${1}\033[0m"
}

# run a command and exit if it fails
function run_cmd {
  echo "RUNNING: $1"
  echo $1
  ec=$?
  if [ $ec -ne 0 ] ; then
    print_bold "ERROR: '$1' returned $?"
    exit $ec
  fi
}

declare -A TOOLS=(
  [torstyverse]="meningotype --version,lissero --version,shovill --version,spades.py -v,skesa --version,mlst --version,prokka --version,snp-dists -v,ngmaster --version,emmtyper --version,csvtk version"
  [seqquality]="seqkit --help,fastp --help,csvtk version"
  [relationships]="kraken2 --version,gubbins -h,mash --version,coresnpfilter --version,iqtree --version,quicktree -v,VeryFastTree --help,gotree version,csvtk version,ska --version"
  [snippy]="snippy --version,csvtk version"
  [mob_suite]="mob_recon --version,csvtk version"
  [panaroo]="panaroo --version,csvtk version"
  [ectyper]="ectyper --version,csvtk version"
  [kleborate]="kleborate --version,csvtk version"
  [stype]="sistr --version,stype --version,csvtk version"
  [tamr]="abritamr --version,tbtamr --version,csvtk version"
  [sonneitype]="mykrobe --version,csvtk version"
  [classify-pangenome]="R --version"
)

# MAIN LOOP OVER ALL ENVS
    
for tool in ${!TOOLS[@]}; do
    print_bold "Setting up: $tool"

    if [[ $FORCE == "true" ]] ; then
        run_cmd "rm -rf $BOHRA_CONDA_ENVS/$tool"
    fi

    if [[ $ACTION == "install" ]]; then
        run_cmd "$INSTALLER env create -p $BOHRA_CONDA_ENVS/$1 -f $ENVSDIR/$1.yml"
    fi

    tests=${TOOLS[$tool]}
    IFS=',' read -r -a cmds <<< "$tests"
    for i in "${!cmds[@]}"; do 
      cmd=${cmds[$i]}
      run_cmd "conda run -p $BOHRA_CONDA_ENVS/$tool $cmd"
      #echo "$tool => $cmd"
    done
done
        
