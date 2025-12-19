#!/usr/bin/env bash

# Safery first!
set -e -u -o pipefail
# set -x

# Issue #133 - Redurectr akk stederr ti stdout 
#exec 2>&1

# check parameters
if [ "$#" -ne 5 ]; then
	EXE=$(basename $0)
	echo "USAGE: $EXE <envsdirr=path> <action=install|check> <force=true|false> <tool=all|toolname> <config=path>"
	exit 255
fi

YAML_DIR=$1
ACTION=$2
FORCE=$3
TOOL=$4
CONFIG=$5

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

declare -A TOOLS
# echo $TOOLS

while IFS= read -r -d '' key && IFS= read -r -d '' value; do
    k=$(echo $key | xargs)
    v=$(echo $value | xargs)
    TOOLS["$k"]="$v"
done < <(jq -r 'to_entries[] | .key + "\u0000" + .value + "\u0000"' "$CONFIG")


# if TOOL is "all", use all tools


# echo "${TOOLS[$TOOL]}"
if [[ "$TOOL" != "all" ]]; then
    
    TMP_TOOLS="${TOOLS[$TOOL]}"
    echo "Selected tools: ${TMP_TOOLS}"
    TOOLS=()
    TOOLS[$TOOL]=$TMP_TOOLS
    
fi

# # MAIN LOOP OVER ALL ENVS
    
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
    # echo "$tool"

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
