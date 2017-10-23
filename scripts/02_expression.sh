#!/bin/bash -l

# Get input parameters
FASTQDIR=$1
EXPRDIR=$2
SAMPLE=$3
LAYOUT=$4
REF=$5

# Specify options based on read layout
if [ "$LAYOUT" == "PAIRED" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
    FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}
    READS="-1$FASTQ1 -2$FASTQ2"

elif [ "$LAYOUT" == "SINGLE" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
    READS="-r$FASTQ1"

fi

# Load modules
module load bioinfo-tools Salmon/0.8.2

# Run Salmon
salmon quant \
    --index $REF \
    --output $EXPRDIR \
    --gcBias \
    --libType A \
    "$READS"
