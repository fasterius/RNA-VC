#!/bin/bash -l

# Get input parameters
WORKDIR=$1
SAMPLE=$2
LAYOUT=$3
GROUP=$4
GEN_REF=$5
STAR_REF=$6

# Get directory for fastq-files and junction files
FASTQDIR=$(echo "$WORKDIR" | sed 's/alignment/fastq/g')
JUNCDIR=$(echo "$WORKDIR" | sed 's/alignment/junctions/g')

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

# Collect junctions
JUNCTIONS=$(find $JUNCDIR/junctions -name "*junctions.tab" | xargs)

# Second pass alignment (with collected first pass junctions across samples)
star --genomeDir $REF_STAR \
    --readFilesIn $FASTQ1 $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN 16 \
    --outSAMattrRGline ID:$SAMPLE LB:$SAMPLE PL:Illumina SM:$SAMPLE \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbFileChrStartEnd $JUNCTIONS \
    --outFileNamePrefix $WORKDIR/pass2/

# Move alignment file
mv $WORKDIR/pass2/Aligned.out.bam > $WORKDIR/${SAMPLE}.bam

# Remove temporary files
rm -r $WORKDIR/pass2

# # Move log file
# cat $WORKDIR/pass2/Log.out > $WORKDIR/log.alignment.pass2.$SAMPLE.txt
