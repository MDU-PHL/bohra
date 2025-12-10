#!/usr/bin/env bash
set -e -u -o pipefail

ENV_PREFIX=$CONDA_PREFIX
GET_DB=$1

# resets to base env
# eval "$(conda shell.bash hook)"
source "$(conda info --base)/etc/profile.d/conda.sh"

DB_VARS=("KRAKEN2_DEFAULT_DB" "BOHRA_PUBMLST_DB" "BOHRA_BLAST_DB" "BOHRA_MOBSUITE_DB")
NONESSENTIAL_VARS=("BOHRA_PUBMLST_DB" "BOHRA_BLAST_DB" "BOHRA_MOBSUITE_DB")
# KRAKEN2_ESSENTIAL_FILES=("hash.k2d" "opts.k2d" "taxo.k2d" "tree.k2d" "names.dmp" "nodes.dmp")
declare -A DB_URL=(
    [1]="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz",
    [2]="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz",
    [3]="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz",
    [4]="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20250402.tar.gz",
    [5]="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20250402.tar.gz",
    [6]="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20250402.tar.gz",
    [7]="https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20250609.tar.gz"

)

if [[ "$GET_DB" =~ "get" ]]; then
    echo "Do you need to download database for kraken2? (y/n)"
    read -r response
    if [[ "$response" == "y" ]]; then
        echo "Please select a database to download:"
        echo "1 - Standard size ~87GB"
        echo "2 - Standard-16 size ~16GB"
        echo "3 - Standard-8 size ~8GB"
        echo "4 - PlusPF size ~94GB"
        echo "5 - PlusPF-16 size ~16GB"
        echo "6 - PlusPF-8 size ~8GB"
        echo "7 - GTDB v226 size ~644GB"
        read -r db_choice
        db_url=${DB_URL[$db_choice]}
        if [[ -z "$db_url" ]]; then
            echo "Invalid choice. Exiting."
            exit 1
        fi
        echo "Downloading database from $db_url"
        wget -q "$db_url" 
        echo "Database downloaded successfully."
        echo "Setting up the database environment variable BOHRA_KRAKEN2_DB"
        echo "Please provide the path where you want to extract the database:"
        read -r db_path
        mkdir -p "$db_path"
        tar -xzf "$(basename "$db_url")" -C "$db_path"
        echo "Database extracted to $db_path"
        conda activate $ENV_PREFIX && conda env config vars set KRAKEN2_DEFAULT_DB="$db_path"
        echo "KRAKEN2_DEFAULT_DB is set to $db_path"
    echo "Reactivating environment to apply changes"
    conda activate $ENV_PREFIX
    else
        echo "Skipping database download."
    fi

fi


for db in "${DB_VARS[@]}"; do
    echo "Checking if $db is set"
        if [ -v "${db}" ]; then
            echo "$db is set and not empty."
        else

            if [[ " ${NONESSENTIAL_VARS[@]} " =~ " $db " ]]; then
                warning_message="It seems that the environment variable $db is not set. If you do not set this variable the pipeline will default to the databases installed in the relevant conda environments or you will need to supply one at the time of running."
            else
                warning_message="It seems that the environment variable $db is not set. You must have at least one species database for the pipeline to be able to detect species from your sequences - if you do not have an environment variable set you will need to provide a database path at the time of running."
            fi
            
            echo "$db is unset or empty. $warning_message"
            if [[ "$GET_DB" =~ "get" ]];then
                echo "Would you like to set it now? (y/n)"
                read -r response
                if [[ "$response" == "y" ]]; then
                    echo "Please enter a value for $db:"
                    read -r value
                    echo "Setting a env var for $db with value: $value"
                    conda activate $ENV_PREFIX && conda env config vars set $db="$value"
                    echo "Reactivating env"
                    conda activate $ENV_PREFIX
                    echo "$db"
                else
                    echo "$db will remain unset."
                fi
            fi
        fi
        echo "Current value of $db: ${!db}"
done

# add in a check to make sure that the kraken2 files are present in the database path and not empty
echo "Finished checking and setting environment variables for databases."

