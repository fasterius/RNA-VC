#!/bin/bash -l

# Get input parameters
WORKDIR=$1
SAMPLE=$2
LAYOUT=$3
GROUP=$4
GEN_REF=$5
STAR_REF=$6

# Get directory for fastq-files
FASTQDIR=$(echo "$WORKDIR" | sed 's/junctions/fastq/g')

# Get read layout and FASTQ input files
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
module load bioinfo-tools star/2.5.3a samtools/1.5

# First pass alignment
star --genomeDir $REF_STAR \
    --readFilesIn $FASTQ1 $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN 16 \
    --outFileNamePrefix $WORKDIR/pass1/ \
    -outSAMmode None

# Move junctions
mv $WORKDIR/pass1/SJ.out.tab $WORKDIR/${SAMPLE}.junctions.tab

# Remove temporary files
rm -r $WORKDIR/pass1

# # Move log file
# cat $WORKDIR/pass1/Log.out > $WORKDIR/log.alignment.pass1.$SAMPLE.txt
