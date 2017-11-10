#!/bin/bash -l

# Get input parameters
FASTQDIR=$1
ALIGNDIR=$2
SAMPLE=$3
LAYOUT=$4
STAR_REF=$5
THREADS=$6
LOGFILE=$7

# Create directory for second pass alignment
TEMP=$ALIGNDIR/temp
mkdir -p $TEMP

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
    --outSAMattrRGline ID:$SAMPLE LB:$SAMPLE PL:Illumina SM:$SAMPLE \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --outFileNamePrefix $TEMP/

# Move alignment file
mv $TEMP/Aligned.sortedByCoord.out.bam $ALIGNDIR/${SAMPLE}.bam

# Append logs
cat $TEMP/Log.out >> $LOGFILE
cat $TEMP/Log.final.out >> $LOGFILE

# Remove temporary files
rm -r $TEMP

# Remove FASTQ-files
rm -r $FASTQDIR
