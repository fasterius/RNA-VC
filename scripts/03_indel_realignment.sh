#!/bin/bash -l

# Bash strict mode
set -euo pipefail
IFS=$'\n\t'

# Input parameters 
ALIGNDIR=$1
REALIGNDIR=$2
SAMPLE=$3
REF=$4
PICARD=$5
GATK=$6
KNOWNINDELS=$7

# Create subdirectories
WORKDIR=$REALIGNDIR/workdir
mkdir -p $REALIGNDIR

# Mark duplicates
java -Xmx5G -jar $PICARD/picard.jar MarkDuplicates \
    INPUT=$ALIGNDIR/${SAMPLE}.bam \
    OUTPUT=$WORKDIR/${SAMPLE}.dedupped.bam  \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    METRICS_FILE=$REALIGNDIR/dedup.metrics.${SAMPLE}.txt

# Remove original bam file
rm $ALIGNDIR/${SAMPLE}.bam

# Index intermediate bam file
samtools index $WORKDIR/${SAMPLE}.dedupped.bam

# Split'N'Trim and reassign mapping qualities
java -Xmx5G -jar $GATK/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R $REF \
    -I $WORKDIR/${SAMPLE}.dedupped.bam \
    -o $WORKDIR/${SAMPLE}.dedupped.split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

# Remove intermediate files
rm $WORKDIR/${SAMPLE}.dedupped.bam*

# Indel realignment
java -Xmx5G -jar $GATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $REF \
    -I $WORKDIR/${SAMPLE}.dedupped.split.bam \
    -known $KNOWNINDELS \
    -o $WORKDIR/realignment_targets.list

java -Xmx5G -jar $GATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $REF \
    -I $WORKDIR/${SAMPLE}.dedupped.split.bam \
    -targetIntervals $WORKDIR/realignment_targets.list \
    -known $KNOWNINDELS \
    -o $REALIGNDIR/${SAMPLE}.bam

# Remove intermediate files
rm -r $WORKDIR

# Index bam file
samtools index $REALIGNDIR/${SAMPLE}.bam
