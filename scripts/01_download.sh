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
    -X 100000 \
    --outdir $FASTQDIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    "$SPLIT" \
    -v \
    $SAMPLE

# Rename paired-end read file names (if applicable)
if [ "$LAYOUT" == "PAIRED" ]; then

    mv $FASTQDIR/${SAMPLE}_1.fastq.gz $FASTQDIR/${SAMPLE}.fastq_1.gz
    mv $FASTQDIR/${SAMPLE}_2.fastq.gz $FASTQDIR/${SAMPLE}.fastq_2.gz
fi

# Delete SRA file (if it exists)
if [ -f $CACHEDIR/${SAMPLE}.sra ]; then
    rm $CACHEDIR/${SAMPLE}.sra
fi

# Delete SRA cache file (if it exists)
if [ -f $CACHEDIR/${SAMPLE}.sra.cache ]; then
    rm $CACHEDIR/${SAMPLE}.sra.cache
fi
