#!/bin/bash -l

# Get input parameters
SAMPLE=$1
LAYOUT=$2
REF=$3

# Working directory
FASTQDIR=data/fastq/$SAMPLE
EXPRDIR=data/expression/$SAMPLE

# Estimate transcript expression with Kallisto
if [ "$LAYOUT" == "PAIRED" ]; then
    FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
    FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}
elif [ "$LAYOUT" == "SINGLE" ]; then
    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
else
    echo "Invalid read layout $LAYOUT; aborting."
    exit 1
fi

# Load modules
module load bioinfo-tools kallisto/0.43.1

# Run Kallisto
kallisto quant \
    -i $REF \
    -t 1 \
    -b 0 \
    -o $EXPRDIR \
    $FASTQ1 $FASTQ2
