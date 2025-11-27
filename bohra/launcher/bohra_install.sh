#!/bin/bash
# Creates conda environments required for bohra

BOHRA_CONDA_ENVS=$CONDA_PREFIX/conda_envs
ENVS_FILES=$1
INSTALL=$2
FORCE_REINSTALL=$3
# abort if any step fails
set -e
# resets to base env
eval "$(conda shell.bash hook)"
echo "Checking your conda version"
condav=$(conda -V | grep -Eo '[0-9]+([.][0-9]+)')
echo Version $condav

# check that mamba is installed
echo Checking which package installer to use
if [ -x "$(which mamba)" ]
then 
    echo "mamba is present and will be used to install dependencies" 
    INSTALLER="mamba"
#   return 0 
else 
    echo "mamba is not installed - will have to use conda - please be patient this may take longer than expected"
    INSTALLER="conda" 
#   return 1 
fi


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

function install_tool(){
    
    echo Running $INSTALLER env create -p $BOHRA_CONDA_ENVS/$1 -f $ENVS_FILES/$1.yml
    time $INSTALLER env create -p $BOHRA_CONDA_ENVS/$1 -f $ENVS_FILES/$1.yml
    echo "Will check $2 setup was successful"
    checkinstalled=$(check_installation $BOHRA_CONDA_ENVS/$1 $1)
    if [[ $checkinstalled -ne 0 ]] 
        then
            echo "There was an error installing $1. Please check the error messages above"
            echo -e "\e[1mExiting now. Please contact the Bohra developers if you need help\e[0m"
            exit 1
    else
        echo "$1 has been set up successfully"
    fi
}

    
for key in ${!TOOLS[@]};do
    if [[ $FORCE_REINSTALL == "true" && $INSTALL == "install" ]]
        then
        echo "Removing existing environment for $key"
        echo "Running rm -rf $BOHRA_CONDA_ENVS/$key"
        rm -rf $BOHRA_CONDA_ENVS/$key
        echo "Will reinstall $key"
        install_tool $key
    else
        echo Checking set up for $key environment
        isinstalled=$(check_installation $BOHRA_CONDA_ENVS/$key $key)
        echo "is installed :$isinstalled"
        if [[ "$isinstalled" != 0 ]]
            then
            echo Looks like $key is not installed. 
            if [[ "$INSTALL" == "install" ]]
                then
                echo Setting up $key environment now
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
        
