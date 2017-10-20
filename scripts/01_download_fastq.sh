#!/bin/bash -l

# Get input parameters
SAMPLE=$1
LAYOUT=$2
REF=$3
CACHEDIR=$4

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
module load bioinfo-tools sratools/2.8.0 kallisto/0.43.1

# Working directory
FASTQDIR=data/fastq/$SAMPLE
EXPRDIR=data/expression/$SAMPLE

# Download FASTQ files
fastq-dump \
    --outdir $FASTQDIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    -v $SAMPLE \
    "$SPLIT" \
        > $FASTQDIR/log.download_fastq.${SAMPLE}.txt 2>&1

# Rename single-end files
if [ "$LAYOUT" == "SINGLE" ]; then
    mv $FASTQDIR/${SAMPLE}.fastq.gz $FASTQDIR/${SAMPLE}_0_fastq.gz
fi

# Delete SRA file (if existing)
if [ -f $CACHEDIR/${SAMPLE}.sra ]; then
    rm $CACHEDIR/${SAMPLE}.sra
fi

# Delete SRA cache file (if existing)
if [ -f $CACHEDIR/${SAMPLE}.sra.cache ]; then
    rm $CACHEDIR/${SAMPLE}.sra.cache
fi

# Estimate transcript expression with Kallisto
if [ "$LAYOUT" == "PAIRED" ]; then
        FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
        FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}
    else
        FASTQ1=$FASTQDIR/${SAMPLE}_0.fastq.gz
        FASTQ2=""
fi

# Run Kallisto
kallisto quant \
    -i $REF \
    -t 1 \
    -b 0 \
    -o $EXPRDIR \
    $FASTQ1 $FASTQ2 \
        > $EXPRDIR/log.kallisto.${SAMPLE}.txt 2>&1
