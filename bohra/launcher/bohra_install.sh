#!/usr/bin/env bash
# Creates conda environments required for bohra
set -e

BOHRA_CONDA_ENVS=$CONDA_PREFIX/conda_envs
ENVS_FILES=$1
INSTALL=$2
FORCE_REINSTALL=$3

INSTALLER="conda"

# resets to base env
eval "$($INSTALLER shell.bash hook)"
$INSTALLER -V

# check that mamba is installed
echo Checking which package installer to use
if [ -x "$(which mamba)" ] ; then
    INSTALLER="mamba"
fi
echo "Using installer: $INSTALLER"

# declare TOOLS=(mash assemblers seqtk seqkit trees prokka snippy mlst kraken2 mob_suite panaroo abritamr gubbins snpdists ectyper emmtyper kleborate lissero meningotype ngmaster stype tbtamr datasmryzr shigapass sonneitype ska2 cluster classify-pangenome fastp)
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
[tamr]="abritamr --version && tbtamr --version && csvtk version"
[sonneitype]="mykrobe --version,csvtk version"
[classify-pangenome]="R --version"
)


function check_installation(){
    # echo "conda activate $1 && ${TOOLS[$2]}"
    IFS=',' read -r -a deps <<< "${TOOLS[$2]}"
    # echo deps array: ${deps[@]}
    found_error=0
    for dep in "${deps[@]}"; do
        # echo "Checking: $dep"
        conda activate $1 && $dep > /dev/null 2>&1
        if [[   "$?" -ne 0 ]] 
            then
                echo -e "\e[1mError: $dep failed\e[0m"
                found_error=1
                # return
        fi
    done
    echo $found_error
    

}

# $1 = the name of the env/tool
function install_tool() {
    icmd="$INSTALLER env create -p $BOHRA_CONDA_ENVS/$1 -f $ENVS_FILES/$1.yml"
    echo "Running: $icmd"
    $icmd
    if [[ $? -ne 0 ]]; then
            echo "FAILED: $icmd"
            exit 2
    else
        echo "$1 has been set up successfully"
    fi
}

    
# MAIN LOOP OVER ALL ENVS
    
for key in ${!TOOLS[@]}; do
    if [[ $FORCE_REINSTALL == "true" && $INSTALL == "install" ]] ; then
        echo "Removing existing environment for $key"
        echo "Running rm -rf $BOHRA_CONDA_ENVS/$key"
        rm -rf $BOHRA_CONDA_ENVS/$key
        echo "Will reinstall $key"
        install_tool $key
    else
        echo "Checking set up for $key environment"
        isinstalled=$(check_installation $BOHRA_CONDA_ENVS/$key $key)
        # echo "is installed :$isinstalled"
        if [[ "$isinstalled" != 0 ]] ; then
            echo "Looks like $key is not installed."
            if [[ "$INSTALL" == "install" ]]; then
                echo "Setting up $key environment now"
                install_tool $key
            else
                echo -e "\e[1m$key is not installed\e[0m"
                echo -e "\e[1mYour installation is incomplete.Exiting now. Please contact the Bohra developers if you need help\e[0m"
                exit 1
            fi
        else
            echo "$key is set up correctly"
        fi
    fi  
done
        
