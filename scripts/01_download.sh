#!/bin/bash -l

# Get input parameters
WORKDIR=$1
SAMPLE=$2
LAYOUT=$3
REF=$4
CACHEDIR=$5

# Get layout from input file name
if [ "$LAYOUT" == "PAIRED" ]; then
    SPLIT="--split-files"
    FASTQTYPE="_1.fastq.gz"
elif [ "$LAYOUT" == "SINGLE" ]; then
    SPLIT=""
    FASTQTYPE=".fastq.gz"
else
    echo "Invalid read layout $LAYOUT; aborting."
    exit 1
fi

# Load modules
module load bioinfo-tools sratools/2.8.0

# Download FASTQ files
fastq-dump \
    --outdir $WORKDIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    -v $SAMPLE \
    "$SPLIT"

# Delete SRA file (if existing)
if [ -f $CACHEDIR/${SAMPLE}.sra ]; then
    rm $CACHEDIR/${SAMPLE}.sra
fi

# Delete SRA cache file (if existing)
if [ -f $CACHEDIR/${SAMPLE}.sra.cache ]; then
    rm $CACHEDIR/${SAMPLE}.sra.cache
fi
