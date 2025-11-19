#!/bin/bash
# Creates conda environments required for bohra

BOHRA_CONDA_ENVS=$CONDA_PREFIX/conda_envs
ENVS_FILES=$1
INSTALL=$2
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

function check_installation(){
    rt=$(conda activate $1 2>&1)
    # echo "$rt"
    if [[ "$rt" =~ .*"Not".*  ]]
    then 
        x=1
    else 
        x=0
    fi
    echo $x
}

declare TOOLS=(mash assemblers seqtk seqkit trees prokka snippy mlst kraken2 mob_suite panaroo abritamr gubbins snpdists ectyper emmtyper kleborate lissero meningotype ngmaster stype tbtamr datasmryzr shigapass sonneitype ska2 cluster classify-pangenome fastp)



if [[ "$INSTALL" == "install" ]]
    then
        for key in ${TOOLS[@]};do
            # echo "$key"
            echo Checking set up for $key
            su=$(check_installation $BOHRA_CONDA_ENVS/$key)
            if [[ $su -eq 1 ]]
                then
                nstll=0
                echo $key can not be found. Now setting up $key
                echo Will now run $INSTALLER env create -p $BOHRA_CONDA_ENVS/$key -f $ENVS_FILES/$key.yml
                $INSTALLER env create -p $BOHRA_CONDA_ENVS/$key -f $ENVS_FILES/$key.yml || nstll=$?
                echo "$nstll"
                echo "Checking installation status for $key"
                su2=$(check_installation $BOHRA_CONDA_ENVS/$key)
                if [[ $su2 -eq 1 ]]
                    then
                    echo "There was an error installing $key. Please check the error messages above"
                    echo -e "\e[1mExiting now. Please contact the Bohra developers if you need help\e[0m"
                    exit 1
                else
                    echo "$key has been successfully installed"
                fi
            fi

        
        done
        
elif [[ "$INSTALL" == "check" ]]
    then
        for key in ${TOOLS[@]};do
            # echo "$key"
            echo Checking set up for $key
            su=$(check_installation $BOHRA_CONDA_ENVS/$key)
            if [[ $su -eq 1 ]]
                then
                echo -e "\e[1m$key is not installed\e[0m"
                echo -e "\e[1mYour installation is incomplete.Exiting now. Please contact the Bohra developers if you need help\e[0m"
                exit 1
            else
                echo "$key is installed"
            fi

            
            done

fi