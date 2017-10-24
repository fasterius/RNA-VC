#!/bin/bash -l

# Get input parameters
FASTQDIR=$1
SAMPLE=$2
LAYOUT=$3
REF=$4
CACHEDIR=$5

# Set paired/single-end reads
if [ "$LAYOUT" == "PAIRED" ]; then

    SPLIT="--split-files"

elif [ "$LAYOUT" == "SINGLE" ]; then

    SPLIT=""
fi

# Download FASTQ files
fastq-dump \
    --outdir $FASTQDIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    "$SPLIT" \
    -v \
    $SAMPLE

# Delete SRA file (if it exists)
if [ -f $CACHEDIR/${SAMPLE}.sra ]; then
    rm $CACHEDIR/${SAMPLE}.sra
fi

# Delete SRA cache file (if it exists)
if [ -f $CACHEDIR/${SAMPLE}.sra.cache ]; then
    rm $CACHEDIR/${SAMPLE}.sra.cache
fi
