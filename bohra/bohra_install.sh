# Creates conda environments required for bohra
# Each environment will be prefaced with 'bohra-' to distinguish them from other
# possible installations that the user has.

ENV_PREFIX="bohra"
# abort if any step fails
set -e
# resets to base env
eval "$(conda shell.bash hook)"
# check that mamba is installed
mamba -V eq 0 
if [ $? -eq 0 ] 
then 
  echo "mamba is installed!! Will no proceed to installations." 
#   return 0 
else 
  echo "mamba is not installed - we will try to install mamba for you" 
  conda install mamba
#   return 1 
fi

function check_installation(){
    rt=$(conda activate $1 2>&1)
    # echo "$rt"
    if [[ "$rt" =~ .*"EnvironmentNameNotFound".* ]]
    then 
        x=1
    else
        x=0
    fi
    echo $x
}

ENVS=($ENV_PREFIX-csvtk $ENV_PREFIX-mash $ENV_PREFIX-spades $ENV_PREFIX-skesa $ENV_PREFIX-seqtk $ENV_PREFIX-seqkit $ENV_PREFIX-iqtree $ENV_PREFIX-quicktree $ENV_PREFIX-prokka $ENV_PREFIX-snippy $ENV_PREFIX-shovill $ENV_PREFIX-mlst $ENV_PREFIX-kraken2 $ENV_PREFIX-mob_suite $ENV_PREFIX-panaroo $ENV_PREFIX-abritamr $ENV_PREFIX-gubbins)
# set up conda envs
echo "Checking conda environments for bohra, if they are not found they will be installed."
# snippy
echo "Checking set up for $ENV_PREFIX-snippy"
su=$(check_installation $ENV_PREFIX-snippy)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-snippy can not be found. Now setting up $ENV_PREFIX-snippy
        mamba create --force -y -n $ENV_PREFIX-snippy snpeff=5.0 snippy=4.4.5 snp-sites
    else
        echo $ENV_PREFIX-snippy is already setup. Nothing left to do
fi
echo "Checking set up for $ENV_PREFIX-snpdists"
su=$(check_installation $ENV_PREFIX-snpdists)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-snpdists can not be found. Now setting up $ENV_PREFIX-snpdists
        mamba create -y -n $ENV_PREFIX-snpdists snp-dists csvtk
    else
        echo $ENV_PREFIX-snpdists is already setup. Nothing left to do
fi
# shovill
echo "Checking set up for $ENV_PREFIX-shovill"
su=$(check_installation $ENV_PREFIX-shovill)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-shovill can not be found. Now setting up $ENV_PREFIX-shovill
        mamba create -y -n $ENV_PREFIX-shovill shovill
    else
        echo $ENV_PREFIX-shovill is already setup. Nothing left to do
fi
# # mlst
echo "Checking set up for $ENV_PREFIX-mlst"
su=$(check_installation $ENV_PREFIX-mlst)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-mlst can not be found. Now setting up $ENV_PREFIX-mlst
        mamba create -y -n $ENV_PREFIX-mlst mlst
    else
        echo $ENV_PREFIX-mlst is already setup. Nothing left to do
fi
# # kraken2
echo "Checking set up for $ENV_PREFIX-kraken2"
su=$(check_installation $ENV_PREFIX-kraken2)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-kraken2 can not be found. Now setting up $ENV_PREFIX-kraken2
        mamba create -y -n $ENV_PREFIX-kraken2 kraken2=2.1.2
    else
        echo $ENV_PREFIX-kraken2 is already setup. Nothing left to do
fi
# mob_suite
echo "Checking set up for $ENV_PREFIX-mob_suite"
su=$(check_installation $ENV_PREFIX-mob_suite)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-mob_suite can not be found. Now setting up $ENV_PREFIX-mob_suite
        mamba create -y -n $ENV_PREFIX-mob_suite mob_suite=3.0.2 numpy=1.21.1
    else
        echo $ENV_PREFIX-mob_suite is already setup. Nothing left to do
        conda activate $ENV_PREFIX-mob_suite
        if [ -d $CONDA_PREFIX/lib/python3.8/site-packages/mob_suite/databases ]
            then
                echo mob_suite database is present
            else
                echo mob_suite databases are not present - will now install at $CONDA_PREFIX/lib/python3.8/site-packages/mob_suite/databases
                mob_init
        fi
fi
# panaroo
echo "Checking set up for $ENV_PREFIX-panaroo"
su=$(check_installation $ENV_PREFIX-panaroo)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-panaroo can not be found. Now setting up $ENV_PREFIX-panaroo
        mamba create -y -n $ENV_PREFIX-panaroo panaroo=1.2.9 csvtk
    else
        echo $ENV_PREFIX-panaroo is already setup. Nothing left to do
fi
# # abritamr
echo "Checking set up for $ENV_PREFIX-abritamr"
su=$(check_installation $ENV_PREFIX-abritamr)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-abritamr can not be found. Now setting up $ENV_PREFIX-abritamr
        mamba create -y -n $ENV_PREFIX-abritamr abritamr=1.0.14
    else
        echo $ENV_PREFIX-abritamr is already setup. Nothing left to do
fi
# # gubbins
# # echo "Setting up $ENV_PREFIX-gubbins"
echo "Checking set up for $ENV_PREFIX-gubbins"
su=$(check_installation $ENV_PREFIX-gubbins)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-gubbins can not be found. Now setting up $ENV_PREFIX-gubbins
        mamba create -y -n $ENV_PREFIX-gubbins gubbins=2.4.1 snp-sites=2.5.1
    else
        echo $ENV_PREFIX-gubbins is already setup. Nothing left to do
fi
# prokka
# echo "Setting up $ENV_PREFIX-prokka"
echo "Checking set up for $ENV_PREFIX-prokka"
su=$(check_installation $ENV_PREFIX-prokka)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-prokka can not be found. Now setting up $ENV_PREFIX-prokka
        mamba create -y -n $ENV_PREFIX-prokka prokka
    else
        echo $ENV_PREFIX-prokka is already setup. Nothing left to do
fi
# quicktree
echo "Checking set up for $ENV_PREFIX-  "
su=$(check_installation $ENV_PREFIX-quicktree)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-quicktree can not be found. Now setting up $ENV_PREFIX-quicktree
        mamba create -y -n $ENV_PREFIX-quicktree quicktree=2.5 newick_utils
    else
        echo $ENV_PREFIX-quicktree is already setup. Nothing left to do
fi
# iqtree
echo "Checking set up for $ENV_PREFIX-iqtree"
su=$(check_installation $ENV_PREFIX-iqtree)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-iqtree can not be found. Now setting up $ENV_PREFIX-iqtree
        mamba create -y -n $ENV_PREFIX-iqtree iqtree=2.1.4 snp-sites=2.5.1
    else
        echo $ENV_PREFIX-iqtree is already setup. Nothing left to do
fi
# seqkit
# echo "Setting up $ENV_PREFIX-seqkit"
echo "Checking set up for $ENV_PREFIX-seqkit"
su=$(check_installation $ENV_PREFIX-seqkit)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-seqkit can not be found. Now setting up $ENV_PREFIX-seqkit
        mamba create -y -n $ENV_PREFIX-seqkit csvtk seqkit=2.1.0
    else
        echo $ENV_PREFIX-seqkit is already setup. Nothing left to do
fi
# seqtk
# echo "Setting up $ENV_PREFIX-seqtk"
echo "Checking set up for $ENV_PREFIX-seqtk"
su=$(check_installation $ENV_PREFIX-seqtk)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-seqtk can not be found. Now setting up $ENV_PREFIX-seqtk
        mamba create -y -n $ENV_PREFIX-seqtk seqtk
    else
        echo $ENV_PREFIX-seqtk is already setup. Nothing left to do
fi
# 
# # skesa
echo "Checking set up for $ENV_PREFIX-skesa"
su=$(check_installation $ENV_PREFIX-skesa)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-skesa can not be found. Now setting up $ENV_PREFIX-skesa
        mamba create -y -n $ENV_PREFIX-skesa skesa
    else
        echo $ENV_PREFIX-skesa is already setup. Nothing left to do
fi
# spades
# echo "Setting up $ENV_PREFIX-spades"
echo "Checking set up for $ENV_PREFIX-spades"
su=$(check_installation $ENV_PREFIX-spades)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-spades can not be found. Now setting up $ENV_PREFIX-spades
        mamba create -y -n $ENV_PREFIX-spades spades=3.15.2
    else
        echo $ENV_PREFIX-spades is already setup. Nothing left to do
fi
# mash
# echo "Setting up $ENV_PREFIX-mash"
echo "Checking set up for $ENV_PREFIX-mash"
su=$(check_installation $ENV_PREFIX-mash)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-mash can not be found. Now setting up $ENV_PREFIX-mash
        mamba create -y -n $ENV_PREFIX-mash mash
    else
        echo $ENV_PREFIX-mash is already setup. Nothing left to do
fi
# csvtk
echo "Checking set up for $ENV_PREFIX-csvtk"
su=$(check_installation $ENV_PREFIX-csvtk)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-csvtk can not be found. Now setting up $ENV_PREFIX-csvtk
        mamba create -y -n $ENV_PREFIX-csvtk csvtk
    else
        echo $ENV_PREFIX-csvtk is already setup. Nothing left to do
fi
# kmc
echo "Checking set up for $ENV_PREFIX-kmc"
su=$(check_installation $ENV_PREFIX-kmc)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-kmc can not be found. Now setting up $ENV_PREFIX-kmc
        mamba create -y -n $ENV_PREFIX-kmc kmc
    else
        echo $ENV_PREFIX-kmc is already setup. Nothing left to do
fi
# stype
echo "Checking set up for $ENV_PREFIX-stype"
su=$(check_installation $ENV_PREFIX-stype)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-stype can not be found. Now setting up $ENV_PREFIX-stype
        mamba create -y -n $ENV_PREFIX-stype sistr_cmd=1.1.1
        conda activate $ENV_PREFIX-stype
        pip3 install git+https://github.com/MDU-PHL/salmonella_typing
    else
        echo $ENV_PREFIX-csvtk is already setup. Nothing left to do
fi
# lissero
echo "Checking set up for $ENV_PREFIX-lissero"
su=$(check_installation $ENV_PREFIX-lissero)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-lissero can not be found. Now setting up $ENV_PREFIX-lissero
        mamba create -y -n $ENV_PREFIX-lissero lissero csvtk
    else
        echo $ENV_PREFIX-lissero is already setup. Nothing left to do
fi
# meningtype
echo "Checking set up for $ENV_PREFIX-meningotype"
su=$(check_installation $ENV_PREFIX-meningotype)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-meningotype can not be found. Now setting up $ENV_PREFIX-meningotype
        mamba create -y -n $ENV_PREFIX-meningotype meningotype csvtk
    else
        echo $ENV_PREFIX-meningotype is already setup. Nothing left to do
fi
# ngmaster
echo "Checking set up for $ENV_PREFIX-ngmaster"
su=$(check_installation $ENV_PREFIX-ngmaster)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-ngmaster can not be found. Now setting up $ENV_PREFIX-ngmaster
        mamba create -y -n $ENV_PREFIX-ngmaster ngmaster csvtk
    else
        echo $ENV_PREFIX-ngmaster is already setup. Nothing left to do
fi
# kleborate
echo "Checking set up for $ENV_PREFIX-kleborate"
su=$(check_installation $ENV_PREFIX-kleborate)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-kleborate can not be found. Now setting up $ENV_PREFIX-kleborate
        mamba create -y -n $ENV_PREFIX-kleborate kleborate csvtk
    else
        echo $ENV_PREFIX-kleborate is already setup. Nothing left to do
fi
# ectyper
echo "Checking set up for $ENV_PREFIX-ectyper"
su=$(check_installation $ENV_PREFIX-ectyper)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-ectyper can not be found. Now setting up $ENV_PREFIX-ectyper
        mamba create -y -n $ENV_PREFIX-ectyper ectyper csvtk
    else
        echo $ENV_PREFIX-ectyper is already setup. Nothing left to do
fi
# emmtyper
echo "Checking set up for $ENV_PREFIX-emmtyper"
su=$(check_installation $ENV_PREFIX-emmtyper)
# echo $su
if [[ $su -eq 1 ]]
    then
        echo $ENV_PREFIX-emmtyper can not be found. Now setting up $ENV_PREFIX-emmtyper
        mamba create -y -n $ENV_PREFIX-emmtyper emmtyper csvtk
    else
        echo $ENV_PREFIX-emmtyper is already setup. Nothing left to do
fi

echo The dependencies for bohra are installed in your default conda path - go forth and analyse!!
echo Please contact us at https://github.com/MDU-PHL/bohra for any issues or concerns