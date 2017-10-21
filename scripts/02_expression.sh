#!/bin/bash -l

# Get input parameters
WORKDIR=$1
SAMPLE=$2
LAYOUT=$3
REF=$4

# Get directory for fastq-files
FASTQDIR=$(echo "$WORKDIR" | sed 's/expression/fastq/g')

echo $FASTQDIR
echo $WORKDIR
echo $SAMPLE
echo $LAYOUT
echo $REF
touch $WORKDIR/${SAMPLE}.abundance.tsv
exit 0

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
