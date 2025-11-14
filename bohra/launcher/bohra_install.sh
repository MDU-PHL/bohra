#!/bin/bash
# Creates conda environments required for bohra
# Each environment will be prefaced with 'bohra-' to distinguish them from other
# possible installations that the user has.

BOHRA_CONDA_ENVS=$CONDA_PREFIX/bohra_envs
ENVS_FILES=$1
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

# declare TOOLS=(mash mash2 seqtk)

# declare -A TOOLS=([$ENV_PREFIX-csvtk]="csvtk=0.33" 
# [$ENV_PREFIX-mash]="mash=2.3 csvtk"
# [$ENV_PREFIX-spades]="spades=3.15.2 csvtk"
# [$ENV_PREFIX-skesa]="skesa csvtk"
# [$ENV_PREFIX-seqtk]="seqtk csvtk"
# [$ENV_PREFIX-seqkit]="seqkit==2.1.0 csvtk"
# [$ENV_PREFIX-iqtree]="iqtree=2.2.0 csvtk gotree=0.4.5"
# [$ENV_PREFIX-quicktree]="quicktree=2.5 gotree=0.4.5 csvtk"
# [$ENV_PREFIX-prokka]="prokka csvtk"
# [$ENV_PREFIX-snippy]="snpeff=5.0 snippy=4.4.5 snp-sites=2.5.1 csvtk"
# [$ENV_PREFIX-shovill]="shovill=1.1.0 csvtk"
# [$ENV_PREFIX-mlst]="mlst=2.19.0 csvtk pandas numpy"
# [$ENV_PREFIX-kraken2]="kraken2=2.1.2 csvtk"
# [$ENV_PREFIX-mob_suite]="mob_suite=3.1.9 csvtk"
# [$ENV_PREFIX-panaroo]="panaroo csvtk"
# [$ENV_PREFIX-abritamr]="abritamr=1.0.19 csvtk"
# [$ENV_PREFIX-gubbins]="csvtk=0.25 gubbins=2.4.1 snp-sites=2.5.1"
# [$ENV_PREFIX-snpdists]="snp-dists=0.8.2 csvtk=0.25"
# [$ENV_PREFIX-ectyper]="ectyper csvtk"
# [$ENV_PREFIX-emmtyper]="emmtyper csvtk"
# [$ENV_PREFIX-kleborate]="kleborate=2.3.2 python=3.9.15 csvtk"
# [$ENV_PREFIX-kmc]="kmc csvtk"
# [$ENV_PREFIX-lissero]="lissero csvtk"
# [$ENV_PREFIX-meningotype]="meningotype csvtk"
# [$ENV_PREFIX-ngmaster]="ngmaster csvtk"
# [$ENV_PREFIX-stype]="sistr_cmd=1.1.1 csvtk"
# [$ENV_PREFIX-tbtamr]="tbtamr=1.0.3 csvtk"
# [$ENV_PREFIX-veryfasttree]="veryfasttree csvtk gotree=0.4.5"
# [$ENV_PREFIX-datasmryzr]="python=3.11 csvtk"
# [$ENV_PREFIX-shigapass]="blast=2.16.0 csvtk"
# [$ENV_PREFIX-sonneitype]="mykrobe pandas csvtk"
# [$ENV_PREFIX-ska2]="ska2 csvtk pandas"
# [$ENV_PREFIX-core-snp-filter]="core-snp-filter csvtk"
# [$ENV_PREFIX-cluster]="numpy pandas scikit-learn csvtk"
# [$ENV_PREFIX-classify-pangenome]="r-data.table=1.17.6 r-ggplot2=3.5.2 r-optparse=1.7.5"
# [$ENV_PREFIX-fastp]="fastp csvtk"
# # [$ENV_PREFIX-tbtamr]="datasmryzr",

# )

DEPS_INSTALLED=0

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
#             if [ "$INSTALL" == "yes" ]
#                 then
#                 echo $key can not be found. Now setting up $key
#                 echo Will now run $INSTALLER create -y -n $key ${TOOLS[$key]}
#                 $INSTALLER create -y -n $key ${TOOLS[$key]}
#                 if [[ "$key" == "$ENV_PREFIX-stype" ]]; then
#                     echo "will install $key from github"
#                     conda activate $key && \
#                     pip3 install 'git+https://github.com/MDU-PHL/salmonella_typing'
#                     # continue
#                 fi
#                 if [[ "$key" == "$ENV_PREFIX-mob_suite" ]]; then
#                     echo "will need to add in db download $key"
#                     echo "Will now initialise the mob_suite databases"
#                     conda activate $key && \
#                     mob_init
#                     # continue
#                 fi

#                 if [[ "$key" == "$ENV_PREFIX-tbtamr" ]]; then
#                     $INSTALLER create --force -y -n $key ${TOOLS[$key]}
#                     echo "Will now finish installing tbtamr"
#                     conda activate $key && \
#                     pip3 install pysam joblib tqdm pydantic requests git+https://github.com/jodyphelan/pathogen-profiler@v4.3.0
#                     # continue
#                 fi
#                 if [[ "$key" == "$ENV_PREFIX-datasmryzr" ]]; then
#                     $INSTALLER create --force -y -n $key ${TOOLS[$key]}
#                     echo "Will now finish installing datasmryzr"
#                     conda activate $key && \
#                     pip3 install git+https://github.com/kristyhoran/datasmryzr
#                     # continue
#                 fi
#             else
#                 echo "$key is not installed. If you want to install it please run bohra check --install-deps. Please note if you do not wish to install, then nextflow will create conda environments for each analsysis. This is not the preferred behaviour."
#                 DEPS_INSTALLED=1
#             fi
#         else
#             echo $key is already setup. Nothing left to do
    fi
# #     # continue
    
    done

# if [[ $DEPS_INSTALLED -eq 1 ]]
#     then
#         echo "One or more dependencies are not installed. Please run bohra check --install-deps to install them"
#         exit 1
#     else
#         echo "All dependencies are installed. Bohra is ready to go!"
#         exit 0
# fi