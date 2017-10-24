#!/bin/bash

# Get input parameters
FASTQDIR=$1
JUNCDIR=$2
SAMPLE=$3
LAYOUT=$4
STAR_REF=$5
THREADS=$6

# Create directory for first pass alignment
PASS1=${FASTQDIR/fastq/alignment}
mkdir -p $PASS1

# Get read layout and FASTQ input files
if [ "$LAYOUT" == "PAIRED" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq_1.gz
    FASTQ2=${FASTQ1/.fastq_1.gz/.fastq_2.gz}

elif [ "$LAYOUT" == "SINGLE" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""

fi

# First pass alignment
star --genomeDir $STAR_REF \
    --readFilesIn $FASTQ1 $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN $THREADS \
    --outFileNamePrefix $PASS1/ \
    --outSAMmode None

# Move junctions
mv $PASS1/SJ.out.tab $JUNCDIR/${SAMPLE}.junctions.tsv

# Remove temporary files
rm -r $PASS1
