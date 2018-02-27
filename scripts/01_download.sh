#!/bin/bash -l

# Get input parameters
FASTQDIR=$1
EXPRDIR=$2
SAMPLE=$3
LAYOUT=$4
SALMON_REF=$5
CACHEDIR=$6

# Set paired/single-end read parameters
if [ "$LAYOUT" == "PAIRED" ]; then

    SPLIT="--split-files"
    FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
    FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}
    READS1="--mates1 $FASTQ1"
    READS2="--mates2 $FASTQ2"

elif [ "$LAYOUT" == "SINGLE" ]; then

    SPLIT=""
    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
    READS1="--unmatedReads $FASTQ1"
    READS2=""
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

# Run Salmon
salmon quant \
    --index $SALMON_REF \
    --output $EXPRDIR \
    --gcBias \
    --libType A \
    $READS1 $READS2

# Rename quantification file
mv $EXPRDIR/quant.sf $EXPRDIR/${SAMPLE}.quant.sf
