#!/usr/bin/env bash

# Safery first!
set -e -u -o pipefail

# Issue #133 - Redurectr akk stederr ti stdout 
exec 2>&1

# check parameters
if [ "$#" -ne 4 ]; then
	EXE=$(basename $0)
	echo "USAGE: $EXE <envsdirr=path> <action=install|check> <force=true|false> [tool=all|toolname]"
	exit 255
fi

YAML_DIR=$1
ACTION=$2
FORCE=$3
TOOL=$4

# echo $TOOL

ENVS_DIR="$CONDA_PREFIX/conda_envs"
INSTALLER="conda"
echo "Using YAML dir: $YAML_DIR"
echo "Testing conda version"
# resets to base env
# eval "$(conda shell.bash hook)"
source "$(conda info --base)/etc/profile.d/conda.sh"

$INSTALLER -V
mkdir -p "$ENVS_DIR"

# check that mamba is installed
if [ -x "$(which mamba)" ] ; then
    INSTALLER="mamba"
fi
echo "Using installer: $INSTALLER"

# print a bold string
print_bold () {
    echo -e "\033[1m${1}\033[0m"
}

disk_space () {
  print_bold "DISK SPACE"
  space=$(df -h .)
  print_bold "$space"
}

# run a command and exit if it fails
run_cmd () {
  local cmd=$1
  echo "RUNNING: $cmd"
  eval "$cmd"
  ec=$?
  if [ $ec -ne 0 ] ; then
    print_bold "ERROR: '$cmd' returned $?"
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

# if TOOL is "all", use all tools

echo "Tools to process: ${TOOL}"
# echo "${TOOLS[$TOOL]}"
if [[ "$TOOL" != "all" ]]; then
    
    TMP_TOOLS="${TOOLS[$TOOL]}"
    echo "Selected tools: ${TMP_TOOLS}"
    TOOLS=()
    TOOLS[$TOOL]=$TMP_TOOLS
    
fi

# MAIN LOOP OVER ALL ENVS
    
for tool in ${!TOOLS[@]}; do
    print_bold "${ACTION}ing : $tool"

    envdir="$ENVS_DIR/$tool"

    if [[ $FORCE == "true" ]] ; then
        run_cmd "rm -rf $envdir"
    fi

    if [[ $ACTION == "install" && ! -d "$envdir" ]]; then
        run_cmd "$INSTALLER env create -p $envdir -f $YAML_DIR/$tool.yml"      
        disk_space
        run_cmd "$INSTALLER clean --all --yes"
        disk_space
    fi

    tests=${TOOLS[$tool]}
    IFS=',' read -r -a cmds <<< "$tests"
    for i in "${!cmds[@]}"; do 
      cmd=${cmds[$i]}
      print_bold "$tool :: $cmd"
      run_cmd "conda run -p $envdir $cmd"
    done
done
        
disk_space
