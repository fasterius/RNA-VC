#!/bin/bash -l
#SBATCH --account b2014056
#SBATCH --partition core
#SBATCH --ntasks 1
#SBATCH --job-name rna_variant_calling
#SBATCH --time 5-00:00:00
#SBATCH --mail-user erikfas@kth.se
#SBATCH --mail-type=FAIL
#SBATCH --error log.rna_variant_calling.err
#SBATCH --output log.rna_variant_calling.out

# Script for variant calling on RNA-seq aligned data
# N.B. Scripts must be run from <PROJECT> main folder with <SAMPLE> as first argument

# Directory structure:
# <PROJECT> /scripts
#           /data    /<SAMPLE>   /00_fastq
#                                /01_expression
#                                /02_alignment
#                                /03_variant.calling
#                                /04_authentication

# BASH strict
set -euo pipefail
IFS=$'\n\t'

# gVCF mode
GVCF="yes"

# Check for missing input
if [ -z ${1+x} ]; then
    echo "ERROR: missing input values; aborting."
    echo "Please provide input SAMPLE as first argument to this script!"
    exit 1
fi

# Sample and log directories
SRR=$1
GSE=$2
ALIGNMENT=data/$GSE/$SRR/02_alignment
WORKDIR=data/$GSE/$SRR/03_variant_calling
REALIGN=$WORKDIR/indel_realignment
BQSR=$WORKDIR/bqsr
CALLS=$WORKDIR/calls
mkdir -p $REALIGN
mkdir -p $BQSR
mkdir -p $CALLS

# Modules
module load bioinfo-tools picard/2.10.3 GATK/3.8-0 samtools snpEff/4.2 bamtools/2.3.0

# Paths
REF=/sw/data/uppnex/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
PICARD=/sw/apps/bioinfo/picard/2.0.1/milou/
GATK=/sw/apps/bioinfo/GATK/3.7/
SNPEFF=/sw/apps/bioinfo/snpEff/4.2/milou/
KNOWNSNPS=/sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/hg38bundle/dbsnp_144.hg38.vcf.gz
KNOWNINDELS=/sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Mark duplicates
java -Xmx7G -jar $PICARD/picard.jar MarkDuplicates \
	INPUT=$ALIGNMENT/${SRR}.bam \
	OUTPUT=$REALIGN/${SRR}.dedupped.bam  \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	METRICS_FILE=$REALIGN/dedup.metrics.${SRR}.txt \
		> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1

# Remove intermediate files
rm $ALIGNMENT/*.bam

# Index intermediate bam file
samtools index $REALIGN/${SRR}.dedupped.bam \
    >> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1

# Split'N'Trim and reassign mapping qualities
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SplitNCigarReads \
	-R $REF \
	-I $REALIGN/${SRR}.dedupped.bam \
	-o $REALIGN/${SRR}.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS \
    	>> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1

# Remove intermediate files
rm $REALIGN/${SRR}.dedupped.bam*

    # Indel realignment
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R $REF \
	-I $REALIGN/${SRR}.dedupped.split.bam \
	-known $KNOWNINDELS \
	-o $REALIGN/realignment_targets.list \
		>> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $REF \
	-I $REALIGN/${SRR}.dedupped.split.bam \
	-targetIntervals $REALIGN/realignment_targets.list \
	-known $KNOWNINDELS \
	-o $BQSR/${SRR}.realigned_indels.bam \
    	>> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1

# Remove intermediate files    
rm $REALIGN/${SRR}.dedupped.split.bam

	# Index bam file
samtools index $BQSR/${SRR}.realigned_indels.bam \
    >> $WORKDIR/log.indel_realignment.${SRR}.txt 2>&1


# Remove intermediate files
rm -r $REALIGN

# First pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R $REF \
	-I $BQSR/${SRR}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
	-knownSites $KNOWNINDELS \
	-o $BQSR/${SRR}.recal.data.table.txt \
		> $WORKDIR/log.bqsr.${SRR}.txt 2>&1

# Second pass covariation modelling (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
    -R $REF \
	-I $BQSR/${SRR}.realigned_indels.bam \
    -knownSites $KNOWNSNPS \
	-knownSites $KNOWNINDELS \
	-BQSR $BQSR/${SRR}.recal.data.table.txt \
	-o $BQSR/${SRR}.post.recal.data.table.txt \
		>> $WORKDIR/log.bqsr.${SRR}.txt 2>&1

# Apply recalibration (BQSR)
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R $REF \
	-I $BQSR/${SRR}.realigned_indels.bam \
    -BQSR $BQSR/${SRR}.recal.data.table.txt \
	-o $BQSR/${SRR}.recal.reads.bam \
		>> $WORKDIR/log.bqsr.${SRR}.txt 2>&1

# Remove intermediate files
rm $BQSR/${SRR}.realigned_indels.bam*

# Variant calling
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R $REF \
	-I $BQSR/${SRR}.recal.reads.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o $CALLS/${SRR}.all.calls.vcf \
	-out_mode EMIT_ALL_CONFIDENT_SITES \
		> $WORKDIR/log.variant_calling.${SRR}.txt 2>&1

# gVCF mode (if applicable)
if [ "$GVCF" == "yes" ]; then
    java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R $REF \
        -I $BQSR/${SRR}.recal.reads.bam \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -emitRefConfidence GVCF \
        -o $WORKDIR/${SRR}.gvcf \
            > $WORKDIR/log.variant_calling_gvcf.${SRR}.txt 2>&1
fi

# Remove intermediate files
rm -r $BQSR

# Separate variants and non-variants for VariantFiltration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REF \
	--variant $CALLS/${SRR}.all.calls.vcf \
	--selectTypeToExclude NO_VARIATION \
	-o $CALLS/${SRR}.variants.vcf \
		>> $WORKDIR/log.variant_calling.${SRR}.txt 2>&1

java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REF \
	--variant $CALLS/${SRR}.all.calls.vcf \
	--selectTypeToInclude NO_VARIATION \
	-o $CALLS/${SRR}.nonvariants.vcf \
        >> $WORKDIR/log.variant_calling.${SRR}.txt 2>&1

# Variant filtration
java -Xmx7G -jar $GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF \
	-V $CALLS/${SRR}.variants.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o $CALLS/${SRR}.variants.filtered.vcf \
        >> $WORKDIR/log.variant_calling.${SRR}.txt 2>&1

# Remove intermediate files
rm $CALLS/${SRR}.variants.vcf

# Concatenate non-variants and filtered variants
java -Xmx7G -cp $GATK/GenomeAnalysisTK.jar \
	org.broadinstitute.gatk.tools.CatVariants \
	-R $REF \
	-V $CALLS/${SRR}.variants.filtered.vcf \
	-V $CALLS/${SRR}.nonvariants.vcf \
	-assumeSorted \
	-out $CALLS/${SRR}.all.calls.filtered.vcf \
        >> $WORKDIR/log.variant_calling.${SRR}.txt 2>&1

# Remove intermediate files
rm $CALLS/${SRR}.variants.filtered.vcf $CALLS/${SRR}.nonvariants.vcf

# Variant annotation
java -Xmx7G -jar $SNPEFF/snpEff.jar GRCh38.82 -stats $CALLS/snpEff.stats.${SRR}.html $CALLS/${SRR}.all.calls.filtered.vcf \
	> $CALLS/${SRR}.all.calls.filtered.annotated.vcf

java -Xmx7G -jar $SNPEFF/SnpSift.jar annotate -id $KNOWNSNPS $CALLS/${SRR}.all.calls.filtered.annotated.vcf \
	> $WORKDIR/${SRR}.vcf

# Remove intermediate files
rm -r $CALLS

