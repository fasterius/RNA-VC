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
    READS1="--mates1 $FASTQ1"
    READS2="--mates2 $FASTQ2"

elif [ "$LAYOUT" == "SINGLE" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
    READS1="--unmatedReads $FASTQ1"
    READS2=""

fi

# Run Salmon
salmon quant \
    --index $REF \
    --output $EXPRDIR \
    --gcBias \
    --libType A \
    $READS1 $READS2

# Rename quantification file
mv $EXPRDIR/quant.sf $EXPRDIR/${SAMPLE}.quant.sf
