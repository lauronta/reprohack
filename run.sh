#!/bin/sh


# check for requirements

# snakemake

hash snakemake 2>/dev/null

if [ $? -ne 0 ]
then
echo 'snakemake : not found'
echo 'make sure you have it in your $PATH or try installing it using apt install snakemake'
exit 1
fi

#singularity

hash singularity 2>/dev/null

if [ $? -ne 0 ]
then
echo 'singularity : not found'
echo 'make sure you have it in your $PATH or try installing it using apt install singularity'
exit 1
fi


echo "this is the base directory : $BASEDIR"
echo "this is the running directory : $(pwd)"


# move to basedir

cd $(dirname $0)
snakemake --use-singularity