#!/bin/bash -l

# Input parameters 
ALIGNDIR=$1
VARIANTDIR=$2
SAMPLE=$3
REF=$4
PICARD=$5
GATK=$6
KNOWNSNPS=$7
KNOWNINDELS=$8
SNPEFF=$9
SNPEFFASSEMBLY=${10}

# Create subdirectories
REALIGN=$VARIANTDIR/indel_realignment
BQSR=$VARIANTDIR/bqsr
CALLS=$VARIANTDIR/variant_calls
mkdir -p $REALIGN $BQSR $CALLS

# Mark duplicates
java -Xmx7G -jar $PICARD/picard.jar MarkDuplicates \
    INPUT=$ALIGNDIR/${SAMPLE}.bam \
    OUTPUT=$REALIGN/${SAMPLE}.dedupped.bam  \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    METRICS_FILE=$REALIGN/dedup.metrics.${SAMPLE}.txt

# Remove original bam file
rm -r $ALIGNDIR

# Index intermediate bam file
samtools index $REALIGN/${SAMPLE}.dedupped.bam

# Split'N'Trim and reassign mapping qualities
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R $REF \
    -I $REALIGN/${SAMPLE}.dedupped.bam \
    -o $REALIGN/${SAMPLE}.dedupped.split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

# Remove intermediate files
rm $REALIGN/${SAMPLE}.dedupped.bam*

# Indel realignment
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $REF \
    -I $REALIGN/${SAMPLE}.dedupped.split.bam \
    -known $KNOWNINDELS \
    -o $REALIGN/realignment_targets.list

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $REF \
    -I $REALIGN/${SAMPLE}.dedupped.split.bam \
    -targetIntervals $REALIGN/realignment_targets.list \
    -known $KNOWNINDELS \
    -o $BQSR/${SAMPLE}.realigned_indels.bam

# Remove intermediate files
rm $REALIGN/${SAMPLE}.dedupped.split.bam

# Index bam file
samtools index $BQSR/${SAMPLE}.realigned_indels.bam

# Remove intermediate files
rm -r $REALIGN

# First pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $REF \
    -I $BQSR/${SAMPLE}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
    -knownSites $KNOWNINDELS \
    -o $BQSR/${SAMPLE}.recal.data.table.txt

# Second pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $REF \
    -I $BQSR/${SAMPLE}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
    -knownSites $KNOWNINDELS \
    -BQSR $BQSR/${SAMPLE}.recal.data.table.txt \
    -o $BQSR/${SAMPLE}.post.recal.data.table.txt

# Apply recalibration (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R $REF \
    -I $BQSR/${SAMPLE}.realigned_indels.bam \
    -BQSR $BQSR/${SAMPLE}.recal.data.table.txt \
    -o $BQSR/${SAMPLE}.recal.reads.bam

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
    -out_mode EMIT_ALL_CONFIDENT_SITES

# Remove intermediate files
rm -r $BQSR

# Separate variants and non-variants for VariantFiltration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $REF \
    --variant $CALLS/${SAMPLE}.all.calls.vcf \
    --selectTypeToExclude NO_VARIATION \
    -o $CALLS/${SAMPLE}.variants.vcf

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $REF \
    --variant $CALLS/${SAMPLE}.all.calls.vcf \
    --selectTypeToInclude NO_VARIATION \
    -o $CALLS/${SAMPLE}.nonvariants.vcf

# Variant filtration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REF \
    -V $CALLS/${SAMPLE}.variants.vcf \
    -window 35 \
    -cluster 3 \
    -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $CALLS/${SAMPLE}.variants.filtered.vcf

# Remove intermediate files
rm $CALLS/${SAMPLE}.variants.vcf

# Concatenate non-variants and filtered variants
java -Xmx7G -cp $GATK/GenomeAnalysisTK.jar \
    org.broadinstitute.gatk.tools.CatVariants \
    -R $REF \
    -V $CALLS/${SAMPLE}.variants.filtered.vcf \
    -V $CALLS/${SAMPLE}.nonvariants.vcf \
    -assumeSorted \
    -out $CALLS/${SAMPLE}.all.calls.filtered.vcf

# Remove intermediate files
rm $CALLS/${SAMPLE}.variants.filtered.vcf $CALLS/${SAMPLE}.nonvariants.vcf

# Variant annotation
java -Xmx7G -jar $SNPEFF/snpEff.jar $SNPEFFASSEMBLY \
    -stats $CALLS/snpEff.stats.${SAMPLE}.html \
    $CALLS/${SAMPLE}.all.calls.filtered.vcf \
	    > $CALLS/${SAMPLE}.all.calls.filtered.annotated.vcf

java -Xmx7G -jar $SNPEFF/SnpSift.jar annotate -id $KNOWNSNPS \
    $CALLS/${SAMPLE}.all.calls.filtered.annotated.vcf \
    	> $VARIANTDIR/${SAMPLE}.vcf

# Compress VCF file
gzip $VARIANTDIR/${SAMPLE}.vcf

# Remove intermediate files
rm -r $CALLS