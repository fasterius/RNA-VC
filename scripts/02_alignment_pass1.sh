#!/bin/bash -l

# Get input parameters
SAMPLE=$1
LAYOUT=$2
GROUP=$3
GEN_REF=$4
STAR_REF=$5

# Directories
FASTQDIR=data/fastq/$SAMPLE
WORKDIR=data/alignment/$SAMPLE
JUNCDIR=data/alignment/junctions/$GROUP
mkdir -p $WORKDIR $JUNCDIR

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

# Set paths
# REF_STAR=/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/
# REF=/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

# First pass alignment
star --genomeDir $REF_STAR \
    --readFilesIn $FASTQ $FASTQ2 \
    --readFilesCommand zcat \
    --runThreadN 16 \
    --outFileNamePrefix $WORKDIR/pass1/

# Move junctions
mv $WORKDIR/pass1/SJ.out.tab $JUNCDIR/${SAMPLE}.junctions.tab

# Move log file
cat $WORKDIR/pass1/Log.out > $WORKDIR/log.alignment.pass1.$SAMPLE.txt
