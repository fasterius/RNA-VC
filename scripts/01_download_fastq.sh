#!/bin/bash -l

# Get info from Snakemake
INFO=$1
SRR=$(echo $INFO | cut -f '/' -f 3)

# Get layout from input file name
LAYOUT=$(echo $INFO | cut -d '/' -f 2)
if [ "$LAYOUT" == "fastq_pe" ]; then
    SPLIT="--split-files"
    FASTQTYPE="_1.fastq.gz"
elif [ "$LAYOUT" == "fastq_se" ]; then
    SPLIT=""
    FASTQTYPE=".fastq.gz"
else
    echo "Invalid read layout $LAYOUT; aborting."
    exit 1
fi

# Load modules
module load bioinfo-tools sratools/2.8.0 kallisto/0.43.1

# Working directory
DATADIR=data/fastq/$SRR
WORKDIR=analysis/expression/$SRR

# Download FASTQ files
fastq-dump \
    --outdir $DATADIR \
    --gzip \
    --skip-technical \
    --readids \
    --clip \
    -v $SRR \
    "$SPLIT" \
        > $DATADIR/log.download_fastq.${SRR}.txt 2>&1

# Delete SRA file (if existing)
if [ -f /proj/b2014056/ncbi/sra/${SRR}.sra ]; then
    rm /proj/b2014056/ncbi/sra/${SRR}.sra
fi

# Delete SRA cache file (if existing)
if [ -f /proj/b2014056/ncbi/sra/${SRR}.sra.cache ]; then
    rm /proj/b2014056/ncbi/sra/${SRR}.sra.cache
fi

# Estimate transcript expression with Kallisto
for FASTQ in $DATADIR/*$FASTQTYPE
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
        -o $WORKDIR \
        $FASTQ $FASTQ2 \
            > $WORKDIR/log.kallisto.${SRR}.txt 2>&1
done
