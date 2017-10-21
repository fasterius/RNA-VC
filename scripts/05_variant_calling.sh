#!/bin/bash -l

# gVCF mode
GVCF="no"

# Check for missing input
if [ -z ${1+x} ]; then
    echo "ERROR: missing input values; aborting."
    echo "Please provide input SAMPLE as first argument to this script!"
    exit 1
fi

# Input parameters 
SAMPLE=$1
REF=$2
PICARD=$3
GATK=$4
KNOWNSNPS=$5
KNOWNINDELS=$6

# Directories
ALIGNMENT=data/alignment/$SAMPLE
WORKDIR=data/variants/$SAMPLE
REALIGN=$WORKDIR/indel_realignment
BQSR=$WORKDIR/bqsr
CALLS=$WORKDIR/calls
mkdir -p $REALIGN
mkdir -p $BQSR
mkdir -p $CALLS

# Load modules
module load bioinfo-tools \
    picard/2.10.3 \
    GATK/3.8-0 \
    samtools/1.5 \
    snpEff/4.2 \
    bamtools/2.3.0

# Mark duplicates
java -Xmx7G -jar $PICARD/picard.jar MarkDuplicates \
	INPUT=$ALIGNMENT/${SAMPLE}.bam \
	OUTPUT=$REALIGN/${SAMPLE}.dedupped.bam  \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	METRICS_FILE=$REALIGN/dedup.metrics.${SAMPLE}.txt \
		> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

# # Remove intermediate files
# rm $ALIGNMENT/*.bam

# Index intermediate bam file
samtools index $REALIGN/${SAMPLE}.dedupped.bam \
    >> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

# Split'N'Trim and reassign mapping qualities
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SplitNCigarReads \
	-R $REF \
	-I $REALIGN/${SAMPLE}.dedupped.bam \
	-o $REALIGN/${SAMPLE}.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS \
    	>> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

# Remove intermediate files
rm $REALIGN/${SAMPLE}.dedupped.bam*

    # Indel realignment
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R $REF \
	-I $REALIGN/${SAMPLE}.dedupped.split.bam \
	-known $KNOWNINDELS \
	-o $REALIGN/realignment_targets.list \
		>> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $REF \
	-I $REALIGN/${SAMPLE}.dedupped.split.bam \
	-targetIntervals $REALIGN/realignment_targets.list \
	-known $KNOWNINDELS \
	-o $BQSR/${SAMPLE}.realigned_indels.bam \
    	>> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

# Remove intermediate files    
rm $REALIGN/${SAMPLE}.dedupped.split.bam

	# Index bam file
samtools index $BQSR/${SAMPLE}.realigned_indels.bam \
    >> $WORKDIR/log.indel_realignment.${SAMPLE}.txt 2>&1

# Remove intermediate files
rm -r $REALIGN

# First pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R $REF \
	-I $BQSR/${SAMPLE}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
	-knownSites $KNOWNINDELS \
	-o $BQSR/${SAMPLE}.recal.data.table.txt \
		> $WORKDIR/log.bqsr.${SAMPLE}.txt 2>&1

# Second pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
    -R $REF \
	-I $BQSR/${SAMPLE}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
	-knownSites $KNOWNINDELS \
	-BQSR $BQSR/${SAMPLE}.recal.data.table.txt \
	-o $BQSR/${SAMPLE}.post.recal.data.table.txt \
		>> $WORKDIR/log.bqsr.${SAMPLE}.txt 2>&1

# Apply recalibration (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R $REF \
	-I $BQSR/${SAMPLE}.realigned_indels.bam \
    -BQSR $BQSR/${SAMPLE}.recal.data.table.txt \
	-o $BQSR/${SAMPLE}.recal.reads.bam \
		>> $WORKDIR/log.bqsr.${SAMPLE}.txt 2>&1

# Remove intermediate files
rm $BQSR/${SAMPLE}.realigned_indels.bam*

# Variant calling
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R $REF \
	-I $BQSR/${SAMPLE}.recal.reads.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o $CALLS/${SAMPLE}.all.calls.vcf \
	-out_mode EMIT_ALL_CONFIDENT_SITES \
		> $WORKDIR/log.variant_calling.${SAMPLE}.txt 2>&1

# gVCF mode (if applicable)
if [ "$GVCF" == "yes" ]; then
    java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R $REF \
        -I $BQSR/${SAMPLE}.recal.reads.bam \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -emitRefConfidence GVCF \
        -o $WORKDIR/${SAMPLE}.gvcf \
            > $WORKDIR/log.variant_calling_gvcf.${SAMPLE}.txt 2>&1
fi

# Remove intermediate files
rm -r $BQSR

# Separate variants and non-variants for VariantFiltration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REF \
	--variant $CALLS/${SAMPLE}.all.calls.vcf \
	--selectTypeToExclude NO_VARIATION \
	-o $CALLS/${SAMPLE}.variants.vcf \
		>> $WORKDIR/log.variant_calling.${SAMPLE}.txt 2>&1

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REF \
	--variant $CALLS/${SAMPLE}.all.calls.vcf \
	--selectTypeToInclude NO_VARIATION \
	-o $CALLS/${SAMPLE}.nonvariants.vcf \
        >> $WORKDIR/log.variant_calling.${SAMPLE}.txt 2>&1

# Variant filtration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-V $CALLS/${SAMPLE}.variants.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o $CALLS/${SAMPLE}.variants.filtered.vcf \
        >> $WORKDIR/log.variant_calling.${SAMPLE}.txt 2>&1

# Remove intermediate files
rm $CALLS/${SAMPLE}.variants.vcf

# Concatenate non-variants and filtered variants
java -Xmx7G -cp $GATK/GenomeAnalysisTK.jar \
	org.broadinstitute.gatk.tools.CatVariants \
	-R $REF \
	-V $CALLS/${SAMPLE}.variants.filtered.vcf \
	-V $CALLS/${SAMPLE}.nonvariants.vcf \
	-assumeSorted \
	-out $CALLS/${SAMPLE}.all.calls.filtered.vcf \
        >> $WORKDIR/log.variant_calling.${SAMPLE}.txt 2>&1

# Remove intermediate files
rm $CALLS/${SAMPLE}.variants.filtered.vcf $CALLS/${SAMPLE}.nonvariants.vcf

# Variant annotation
java -Xmx7G -jar $SNPEFF/snpEff.jar $SNPEFFASSEMBLY \
    -stats $CALLS/snpEff.stats.${SAMPLE}.html \
    $CALLS/${SAMPLE}.all.calls.filtered.vcf \
	    > $CALLS/${SAMPLE}.all.calls.filtered.annotated.vcf

java -Xmx7G -jar $SNPEFF/SnpSift.jar annotate -id $KNOWNSNPS \
    $CALLS/${SAMPLE}.all.calls.filtered.annotated.vcf \
    	> $WORKDIR/${SAMPLE}.vcf

# Remove intermediate files
rm -r $CALLS
