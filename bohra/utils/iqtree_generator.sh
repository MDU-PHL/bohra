#!/bin/bash

# http://www.iqtree.org/doc/Command-Reference
# courtesy of Torsten Seemann. Modified here to account for snp-detector file structure
#-fconst
#Specify a list of comma-separated integer numbers. The number of
#entries should be equal to the number of states in the model (e.g.  4 for
#DNA and 20 for protein).  IQ-TREE will then add a number of constant sites
#accordingly.  For example, -fconst 10,20,15,40 will add 10 constant sites
#of all A, 20 constant sites of all C, 15 constant sites of all G and 40
#constant sites of all T into the alignment.

set -e

FASTA=$1
if [ "x$FASTA" == "x" ]; then
  EXE=$(basename $0)
  echo "Usage: $EXE genome.fasta core.aln prefix max_cpus"
  exit 1
fi

RESULT="-fconst "
for L in Aa Cc Gg Tt ; do
  NUM=$(grep -v '>' "$FASTA" | tr -d -c "[$L]" | wc -c)
  if [ "x$NUM" == "x" ]; then
    echo "Error: could not count '$L' chars in $FASTA"
    exit 2
  fi
  RESULT="$RESULT$NUM,"
done

# remove last comma
RESULT=${RESULT: : -1}

# get core
CORE=$2

PRE=$3

if [ "x$4" == "x" ]; then
    CPUS=$(getconf _NPROCESSORS_ONLN)
else
    CPUS=$4
fi


echo "iqtree $RESULT -m GTR+G4 -bb 1000 -alrt 1000 -ntmax "$CPUS" -nt AUTO -st DNA -s $CORE -pre $PRE"