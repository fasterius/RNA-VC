#!/bin/bash -l

# Bash strict mode
set -euo pipefail
IFS=$'\n\t'

# Input parameters 
REALIGNDIR=$1
VARIANTDIR=$2
NAME=$3
REF=$4
GATK=$5
KNOWNSNPS=$6
KNOWNINDELS=$7
SNPEFF=$8
SNPEFFASSEMBLY=$9
JAVEMEM="-Xmx${10}"
MERGE=${11}

# Create subdirectories
BQSR=$VARIANTDIR/$NAME/bqsr
CALLS=$VARIANTDIR/$NAME/variant_calls
mkdir -p $BQSR $CALLS

# Check for groups and merge as applicable
if [ "$MERGE" == "MERGE" ]; then

    # List all BAM files in realignment directory
    for BAM in $REALIGNDIR/*.bam; do
        echo "$BAM" >> $BQSR/bam_list.txt
    done

    # Merge and index indel realignment bam files
    bamtools merge -list $BQSR/bam_list.txt -out $BQSR/$NAME.merged.bam
    samtools index $BQSR/$NAME.merged.bam

    # Set input BAM
    INPUT_BAM=$BQSR/${NAME}.merged.bam

    # Remove BAM list
    rm $BQSR/bam_list.txt
else

    # Set input without merging groups
    INPUT_BAM=$REALIGNDIR/${NAME}.bam
fi

# First pass covariation modelling (BQSR)
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $REF \
    -I $INPUT_BAM \
    -knownSites $KNOWNSNPS \
    -knownSites $KNOWNINDELS \
    -o $BQSR/${NAME}.recal.data.table.txt

# Second pass covariation modelling (BQSR)
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $REF \
    -I $INPUT_BAM \
    -knownSites $KNOWNSNPS \
    -knownSites $KNOWNINDELS \
    -BQSR $BQSR/${NAME}.recal.data.table.txt \
    -o $BQSR/${NAME}.post.recal.data.table.txt

# Apply recalibration (BQSR)
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R $REF \
    -I $INPUT_BAM \
    -BQSR $BQSR/${NAME}.recal.data.table.txt \
    -o $BQSR/${NAME}.recal.reads.bam

# Remove intermediate files
rm $REALIGNDIR/*.bam

# Variant calling
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $REF \
    -I $BQSR/${NAME}.recal.reads.bam \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
    -o $CALLS/${NAME}.all.calls.vcf \
    -out_mode EMIT_ALL_CONFIDENT_SITES

# Remove intermediate files
rm -r $BQSR

# Separate variants and non-variants for VariantFiltration
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $REF \
    --variant $CALLS/${NAME}.all.calls.vcf \
    --selectTypeToExclude NO_VARIATION \
    -o $CALLS/${NAME}.variants.vcf

java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $REF \
    --variant $CALLS/${NAME}.all.calls.vcf \
    --selectTypeToInclude NO_VARIATION \
    -o $CALLS/${NAME}.nonvariants.vcf

# Variant filtration
java $JAVAMEM -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REF \
    -V $CALLS/${NAME}.variants.vcf \
    -window 35 \
    -cluster 3 \
    -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $CALLS/${NAME}.variants.filtered.vcf

# Remove intermediate files
rm $CALLS/${NAME}.variants.vcf

# Concatenate non-variants and filtered variants
java $JAVAMEM -cp $GATK/GenomeAnalysisTK.jar \
    org.broadinstitute.gatk.tools.CatVariants \
    -R $REF \
    -V $CALLS/${NAME}.variants.filtered.vcf \
    -V $CALLS/${NAME}.nonvariants.vcf \
    -assumeSorted \
    -out $CALLS/${NAME}.all.calls.filtered.vcf

# Remove intermediate files
rm $CALLS/${NAME}.variants.filtered.vcf
rm $CALLS/${NAME}.nonvariants.vcf

# Change chrM to chrMT (if applicable)
cat $CALLS/${NAME}.all.calls.filtered.vcf | sed 's/^chrM/chrMT/' \
    > $CALLS/${NAME}.renamed.vcf

# Variant annotation (snpEff)
java $JAVAMEM -jar $SNPEFF/snpEff.jar $SNPEFFASSEMBLY \
    -stats $CALLS/snpEff.stats.${NAME}.html \
    $CALLS/${NAME}.renamed.vcf \
	    > $CALLS/${NAME}.all.calls.filtered.annotated.vcf

# Variant annotation (SnpSift)
java $JAVAMEM -jar $SNPEFF/SnpSift.jar annotate -id $KNOWNSNPS \
    $CALLS/${NAME}.all.calls.filtered.annotated.vcf \
    	> $VARIANTDIR/${NAME}.vcf

# Compress VCF file
gzip $VARIANTDIR/${NAME}.vcf

# Remove intermediate files
rm -r $VARIANTDIR/$NAME
