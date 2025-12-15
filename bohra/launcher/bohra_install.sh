#!/usr/bin/env bash

# Safery first!
set -e -u -o pipefail

# Issue #133 - Redurectr akk stederr ti stdout 
#exec 2>&1

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

#echo "Using YAML dir: $YAML_DIR"
#echo "Testing conda version"

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
  df -h
}

# run a command and exit if it fails
run_cmd () {
  local cmd=$1
  echo "Running: $cmd"
  eval "$cmd"
  ec=$?
  if [ $ec -ne 0 ] ; then
    print_bold "ERROR: '$cmd' returned $?"
    exit $ec
  fi
  # disk_space
}

# disk_space

declare -A TOOLS=(
  [torstyverse]="meningotype --version,lissero --version,shovill --version,spades.py -v,skesa --version,mlst --version,prokka --version,snp-dists -v,ngmaster --version,emmtyper --version,csvtk version"
  [seqquality]="seqkit version,fastp --version,csvtk version"
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
        run_cmd "$INSTALLER env create --quiet -p $envdir -f $YAML_DIR/$tool.yml"      
        # clean up conda
        #run_cmd "$INSTALLER clean --all --yes --quiet"
        # remove crap
        run_cmd "find $envdir -name '*.pyc' -delete"
        run_cmd "find $envdir -type d -name PyQt5 | xargs rm -fr"
        run_cmd "find $envdir -name ncbi_plasmid_full_seqs.fas -print -delete"
        run_cmd "find $envdir -name src.zip -print -delete"
        run_cmd "find $envdir -name '*.a' -print -delete"
        run_cmd "rm -fr $envdir/share/EMBOSS"
        run_cmd "rm -fr $envdir/{man,include,docs,doc,legal}"
        run_cmd "rm -fr $envdir/lib/libLLVM*"
        disk_space
    fi

    tests=${TOOLS[$tool]}
    IFS=',' read -r -a cmds <<< "$tests"
    numtests=${#cmds[@]}
    for i in "${!cmds[@]}"; do 
      cmd=${cmds[$i]}
      # c=$i + 1
      print_bold "$tool :: $i/$numtests :: $cmd"
      run_cmd "conda run -p $envdir $cmd"
      
    done
done
        
# disk_space
