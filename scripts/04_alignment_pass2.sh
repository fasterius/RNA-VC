#!/bin/bash -l

# Get input parameters
FASTQDIR=$1
JUNCDIR=$2
ALIGNDIR=$3
SAMPLE=$4
LAYOUT=$5
GROUP=$6
STAR_REF=$7
THREADS=$8

# Create directory for second pass alignment
PASS2=$ALIGNDIR/pass2
mkdir -p $PASS2

# Get read layout and FASTQ input files
if [ "$LAYOUT" == "PAIRED" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}_1.fastq.gz
    FASTQ2=${FASTQ1/_1.fastq.gz/_2.fastq.gz}

elif [ "$LAYOUT" == "SINGLE" ]; then

    FASTQ1=$FASTQDIR/${SAMPLE}.fastq.gz
    FASTQ2=""
fi

# Collect junctions
JUNCTIONS=$(find $JUNCDIR -name "*junctions.tsv" | xargs)

# Second pass alignment
star --genomeDir $STAR_REF \
    --readFilesIn $FASTQ1 $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN $THREADS \
    --outSAMattrRGline ID:$SAMPLE LB:$SAMPLE PL:Illumina SM:$SAMPLE \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbFileChrStartEnd $JUNCTIONS \
    --outFileNamePrefix $PASS2/

# Move alignment file
mv $PASS2/Aligned.sortedByCoord.out.bam $ALIGNDIR/${SAMPLE}.bam.tmp

# Remove temporary files
rm -r $PASS2
