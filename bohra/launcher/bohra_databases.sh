#!/bin/bash

ENV_PREFIX=$1


# abort if any step fails
set -e
# resets to base env
eval "$(conda shell.bash hook)"


DB_VARS=("BOHRA_KRAKEN2_DB" "BOHRA_SYLPH_DB" "BOHRA_PUBMLST_DB" "BOHRA_BLAST_DB" "BOHRA_MOBSUITE_DB")

for db in "${DB_VARS[@]}"; do
    echo "Checking if $db is set"
        if [ -v "${db}" ]; then
            echo "$db is set and not empty."
        else
            echo "$db is unset or empty. It is recommended that you set these variables for optimal performance of the bohra pipeline. Would you like to set it now? (y/n)"
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
        echo "Current value of $db: ${!db}"
done

echo "Finished checking and setting environment variables for databases."