#!/bin/bash -l

# Get input parameters
SRR=$1
LAYOUT=$2
REF=$3
CACHEDIR=$4

# Get layout from input file name
if [ "$LAYOUT" == "PAIRED" ]; then
    SPLIT="--split-files"
    FASTQTYPE="_1.fastq.gz"
elif [ "$LAYOUT" == "SINGLE" ]; then
    SPLIT=""
    FASTQTYPE=".fastq.gz"
else
    echo "Invalid read layout $LAYOUT; aborting."
    exit 1
fi

# Load modules
module load bioinfo-tools sratools/2.8.0 kallisto/0.43.1

# Working directory
FASTQDIR=data/fastq/$SRR
EXPRDIR=data/expression/$SRR

# Download FASTQ files
fastq-dump \
    --outdir $FASTQDIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    -v $SRR \
    "$SPLIT" \
        > $FASTQDIR/log.download_fastq.${SRR}.txt 2>&1

# Delete SRA file (if existing)
if [ -f $CACHEDIR/${SRR}.sra ]; then
    rm $CACHEDIR/${SRR}.sra
fi

# Delete SRA cache file (if existing)
if [ -f $CACHEDIR/${SRR}.sra.cache ]; then
    rm $CACHEDIR/${SRR}.sra.cache
fi

# Estimate transcript expression with Kallisto
for FASTQ in $FASTQDIR/*$FASTQTYPE
do
    # Set variable input files for PE/SE reads
        if [ "$LAYOUT" == "PAIRED" ]; then
            FASTQ2=${FASTQ/_1.fastq.gz/_2.fastq.gz}
        else
            FASTQ2=""
        fi

    # Run Kallisto
    kallisto quant \
        -i $REF \
        -t 1 \
        -b 0 \
        -o $EXPRDIR \
        $FASTQ $FASTQ2 \
            > $EXPRDIR/log.kallisto.${SRR}.txt 2>&1
done
