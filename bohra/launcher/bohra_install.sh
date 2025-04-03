#!/bin/bash
# Creates conda environments required for bohra
# Each environment will be prefaced with 'bohra-' to distinguish them from other
# possible installations that the user has.

ENV_PREFIX=$1
UPDATE=$2
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
    if [[ "$rt" =~ .*"EnvironmentNameNotFound".* || "$rt" =~ .*"Could not find".* ]]
    then 
        x=1
    else 
        x=0
    fi
    echo $x
}



declare -A TOOLS=([$ENV_PREFIX-csvtk]="csvtk=0.33" 
[$ENV_PREFIX-mash]="mash=2.3"
[$ENV_PREFIX-spades]="spades"
[$ENV_PREFIX-skesa]="skesa"
[$ENV_PREFIX-seqtk]="seqtk"
[$ENV_PREFIX-seqkit]="seqkit==2.1.0"
[$ENV_PREFIX-iqtree]="iqtree=2.2.0"
[$ENV_PREFIX-quicktree]="quicktree=2.5 newick_utils"
[$ENV_PREFIX-prokka]="prokka"
[$ENV_PREFIX-snippy]="snpeff=5.0 snippy=4.4.5 snp-sites=2.5.1"
[$ENV_PREFIX-shovill]="shovill=1.1.0"
[$ENV_PREFIX-mlst]="mlst=2.19.0"
[$ENV_PREFIX-kraken2]="kraken2=2.1.2"
[$ENV_PREFIX-mob_suite]="mob_suite=3.1.9"
[$ENV_PREFIX-panaroo]="panaroo=1.2.9 csvtk"
[$ENV_PREFIX-abritamr]="abritamr=1.0.19"
[$ENV_PREFIX-gubbins]="csvtk=0.25 gubbins=2.4.1 snp-sites=2.5.1"
[$ENV_PREFIX-snpdists]="snp-dists=0.8.2 csvtk=0.25"
[$ENV_PREFIX-ectyper]="ectyper csvtk"
[$ENV_PREFIX-emmtyper]="emmtyper csvtk"
[$ENV_PREFIX-kleborate]="kleborate csvtk"
[$ENV_PREFIX-kmc]="kmc"
[$ENV_PREFIX-lissero]="lissero csvtk"
[$ENV_PREFIX-meningotype]="meningotype csvtk"
[$ENV_PREFIX-ngmaster]="ngmaster csvtk"
[$ENV_PREFIX-stype]="sistr_cmd=1.1.1 csvtk"
)

for key in "${!TOOLS[@]}";do
    # do echo "$key - ${TOOLS[$key]}"
    echo Checking set up for $key
    su=$(check_installation $key)
    if [[ $su -eq 1 ]]
        then
            echo $key can not be found. Now setting up $key
            echo Will now run $INSTALLER create --force -y -n $key ${TOOLS[$key]}
            $INSTALLER create --force -y -n $key ${TOOLS[$key]}
            if [[ "$key" == "$ENV_PREFIX-stype" ]]; then
                echo "will install $key from github"
                conda activate $key && \
                pip3 install 'git+https://github.com/MDU-PHL/salmonella_typing'
                # continue
            fi
            if [[ "$key" == "$ENV_PREFIX-mob_suite" ]]; then
                echo "will need to add in db download $key"
                echo "Will now initialise the mob_suite databases"
                conda activate $key && \
                mob_init
                # continue
            fi
        else
            echo $key is already setup. Nothing left to do
        fi
    # continue
    
    done

if [[ -z $KRAKEN2_DEFAULT_DB ]]; then
  
  echo KRAKEN2_DEFAUT_DB is undefined. Would you like to download babykraken database?
  read -p "Download babykraken (10MB) (Y/N)? " download
    if [ "$download" != "${download#[Yy]}" ] ;then 
    
        echo Downloading babykraken now.
        curl -L https://github.com/MDU-PHL/babykraken/blob/master/dist/babykraken.tar.gz?raw=true | tar xz
        echo babykraken has been downloaded. Please use --kraken_db $(pwd)/babykraken to run kraken2
        
    else
        
        echo Okey dokey .. hopefully you remember to provide a kraken path when you run bohra
    fi

   
fi

echo The dependencies for bohra are installed in your default conda path - go forth and analyse!!

echo Please contact us at https://github.com/MDU-PHL/bohra for any issues or concerns