#!/bin/bash -l

# Get info from Snakemake
INFO=$1
SRR=$(echo $INFO | cut -d '/' -f 3)

# Get layout from input file name
LAYOUT=$(echo $INFO | cut -d '/' -f 2)
if [ "$LAYOUT" == "fastq_pe" ]; then
    FASTQTYPE="_1.fastq.gz"
elif [ "$LAYOUT" == "fastq_se" ]; then
    FASTQTYPE=".fastq.gz"
else
    echo "Invalid read layout $LAYOUT; aborting."
    exit 1
fi

# Directories
# FASTQDIR=data/$SRR/$SRR/00_fastq
WORKDIR=analysis/alignment/$SRR
mkdir -p $WORKDIR/junctions

# Load modules
module load bioinfo-tools star/2.5.3a samtools/1.5

# Set paths
REF_STAR=/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/
REF=/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

# First pass (individually)
for FASTQ in $FASTQDIR/*$FASTQTYPE
do
    # Paths and names
    BASENAME=$(basename $FASTQ)
    REPLICATE=${BASENAME/$FASTQTYPE/}

    # Set --readFilesIn variable for PE/SE reads
    if [ "$LAYOUT" == "PAIRED" ]; then
        FASTQ2=${FASTQ/_1.fastq.gz/_2.fastq.gz}
    else
        FASTQ2=""
    fi

    # Alignment
	mkdir -p $WORKDIR/pass1
	star --genomeDir $REF_STAR \
		--readFilesIn $FASTQ $FASTQ2 \
		--readFilesCommand zcat \
		--runThreadN 16 \
        --outFileNamePrefix $WORKDIR/pass1/

    # Move junctions
    mv $WORKDIR/pass1/SJ.out.tab $WORKDIR/junctions/junctions.$REPLICATE.tab

    # Log file
    cat $WORKDIR/pass1/Log.out > $WORKDIR/log.alignment.$REPLICATE.txt

    # Delete temporary files
    rm -r $WORKDIR/pass1
done

# Collect junctions
JUNCTIONS=$(find $WORKDIR/junctions -name "junctions.*.tab" | xargs)

# Second pass (with collected first pass junctions across samples)
for FASTQ in $FASTQDIR/*$FASTQTYPE
do
    # Paths and names
    BASENAME=$(basename $FASTQ)
    REPLICATE=${BASENAME/$FASTQTYPE/}

    # Set --readFilesIn variable for PE/SE reads
    if [ "$LAYOUT" == "PAIRED" ]; then
        FASTQ2=${FASTQ/_1.fastq.gz/_2.fastq.gz}
    else
        FASTQ2=""
    fi
    
    # Alignment
	mkdir $WORKDIR/pass2
	star --genomeDir $REF_STAR \
		--readFilesIn $FASTQ $FASTQ2 \
		--readFilesCommand zcat \
		--runThreadN 16 \
		--outSAMattrRGline ID:$REPLICATE LB:$REPLICATE PL:Illumina SM:$SRR \
		--outSAMtype BAM Unsorted \
        --sjdbFileChrStartEnd $JUNCTIONS \
        --outFileNamePrefix $WORKDIR/pass2/

	# Sort and move alignment file
    samtools sort --threads 16 $WORKDIR/pass2/Aligned.out.bam > $WORKDIR/$REPLICATE.bam
	
	# Append log file
	cat $WORKDIR/pass2/Log.out >> $WORKDIR/log.alignment.$REPLICATE.txt

	# Remove temporary files
	rm -r $WORKDIR/pass2
done

# Remove temporary files
rm $FASTQDIR/*.fastq.gz
rm -r $WORKDIR/junctions
