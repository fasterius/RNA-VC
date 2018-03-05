#!/bin/bash -l

# Set bash strict mode
set -euo pipefail
IFS=$'\n\t'

# Get input parameters
FASTQDIR=$1
ALIGNDIR=$2
SAMPLE=$3
RGSM=$4
LAYOUT=$5
STAR_REF=$6
THREADS=$7
LOGFILE=$8

# Create directory for second pass alignment
WORKDIR=$ALIGNDIR/$SAMPLE
mkdir -p $WORKDIR

# Get read layout and FASTQ input files
if [ "$LAYOUT" == "PAIRED" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
    FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}

elif [ "$LAYOUT" == "SINGLE" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
fi

# Perform alignment
star --genomeDir $STAR_REF \
    --readFilesIn $FASTQ1 $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN $THREADS \
    --outSAMattrRGline ID:$SAMPLE LB:$SAMPLE PL:Illumina SM:$RGSM \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --outFileNamePrefix $WORKDIR/

# Move alignment file
mv $WORKDIR/Aligned.sortedByCoord.out.bam $ALIGNDIR/${SAMPLE}.bam

# Append logs
cat $WORKDIR/Log.out >> $LOGFILE
cat $WORKDIR/Log.final.out >> $LOGFILE

# Remove workdirorary files
rm -r $WORKDIR
rm -r $FASTQDIR
