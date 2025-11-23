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
declare -A TOOLS=([mash]="mash --version && csvtk --version"
[assemblers1]="shovill --version && spades.py -v && skesa --version && csvtk --version"
# [seqtk]="seqtk --version && csvtk --version"
# [seqkit]="csvtk --version" # seqkit does not have a --version flag
# [trees]="coresnpfilter --version && iqtree --version && quicktree -v && VeryFastTree --help && gotree version && csvtk --version"
# [prokka]="prokka --version && csvtk --version"
# [snippy]="snippy --version && csvtk --version"
# [mlst]="mlst --version && csvtk --version"
# [kraken2]="kraken2 --version && csvtk --version"
# [mob_suite]="mob_recon --version && csvtk --version"
# [panaroo]="panaroo --version && csvtk --version"
# [abritamr]="abritamr --version && csvtk --version"
# [gubbins]="gubbins -h && csvtk --version"
# [snpdists]="snp-dists --version && csvtk --version"
# [ectyper]="ectyper --version && csvtk --version"
# [emmtyper]="emmtyper --version && csvtk --version"
# [kleborate]="kleborate --version && csvtk --version"
# [lissero]="lissero --version && csvtk --version"
# [meningotype]="meningotype --version && csvtk --version"
# [ngmaster]="ngmaster --version && csvtk --version"
# [stype]="sistr --version && stype --version && csvtk --version"
# [tbtamr]="tbtamr --version && csvtk --version"
# [datasmryzr]="datasmryzr --help && csvtk --version"
# [shigapass]="blastn -version && csvtk --version"
# [sonneitype]="mykrobe --version && csvtk --version"
# [ska2]="ska  --version && csvtk --version"
# [cluster]="csvtk --version"
# [classify-pangenome]="R --version"
# [fastp]="fastp --version && csvtk --version"
)


function check_installation(){
    rt=$(conda activate $1 && ${TOOLS[$2]} 2>&1 )
    echo "$?"
    
}

function install_tool(){
    echo Installing $1
    time $INSTALLER env create -p $BOHRA_CONDA_ENVS/$1 -f $ENVS_FILES/$1.yml
    
    echo "Will confirm that $1 has been installed"
    checkinstalled=$(check_installation $BOHRA_CONDA_ENVS/$1 $1)
    if [[ $checkinstalled -eq 1 ]] 
        then
            echo "There was an error installing $1. Please check the error messages above"
            echo -e "\e[1mExiting now. Please contact the Bohra developers if you need help\e[0m"
            exit 1
    else
        echo "$1 installed"
    fi
}

    
for key in ${!TOOLS[@]};do

    echo Checking set up for $key
    isinstalled=$(check_installation $BOHRA_CONDA_ENVS/$key $key)
    
    if [[ "$isinstalled" != 0 ]]
        then
        echo Looks like $key is not installed. 
        if [[ "$INSTALL" == "install" ]]
            then
            echo Installing $key now
            install_tool $key
        else
            echo -e "\e[1m$key is not installed\e[0m"
            echo -e "\e[1mYour installation is incomplete.Exiting now. Please contact the Bohra developers if you need help\e[0m"
            exit 1
        fi
    elif [[ "$FORCE_REINSTALL" == "true" && "$isinstalled" == 0 ]]
        then
        echo "Force reinstall is set to true. Will reinstall $key"
        install_tool $key
    else
        echo "$key is installed"
    fi
        
done
        
