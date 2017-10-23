#!/usr/bin/env python3

# Packages
import pandas as pd
import re

# Config file
configfile: "config.yaml"

# Read metadata
metadata_file = config["metadata"]
metadata = pd.read_table(metadata_file, encoding = 'ISO-8859-1')

# Get relevant column names from config
study_col = config["study_col"]
group_col = config["group_col"]
sample_col = config["sample_col"]
layout_col = config["layout_col"]

# Get studies, samples and groups from metadata
STUDIES = metadata[study_col]
STUDIES_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][study_col]
STUDIES_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][study_col]
GROUPS = metadata[group_col]
GROUPS_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][group_col]
GROUPS_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][group_col]
SAMPLES = metadata[sample_col]
SAMPLES_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][sample_col]
SAMPLES_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][sample_col]

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Strings for base directories
base = 'data/{study}/{group}/'
tmpdir = '.tmp_snake'
tmp = tmpdir + '/' + base

# Rule: final outputs
rule all:
    input:
        expand(base + 'expression/{sample}/{sample}.quant.sf', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES),
        expand(base + 'variants/{sample}.vcf', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES)
    shell:
        'rm -rf {tmpdir}'

# Rule: clean
rule clean:
    shell:
        'rm -rf data {tmpdir}'

# Rule: download raw data
rule download_SE:
    output:
        # tmp + 'fastq/{sample}/{sample}.fastq.gz'
        expand(tmp + 'fastq/{sample}/{sample}.fastq.gz', zip,
            study = STUDIES_SE, group = GROUPS_SE, sample = SAMPLES_SE),
        # expand(tmp + 'fastq/{sample}/{sample}_1.fastq.gz', zip,
            # study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE),
        # expand(tmp + 'fastq/{sample}/{sample}_2.fastq.gz', zip,
            # study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE)
    # params:
        # layout = lambda wildcards:
            # get_metadata(wildcards.sample, layout_col),
    log:
        # base + 'logs/01_download.{sample}.log'
        expand(base + 'logs/01_download.{sample}.log', zip,
            study = STUDIES_SE, group = GROUPS_SE, sample = SAMPLES_SE)
    shell:
        """
        FASTQ=$(basename {output})
        SAMPLE=$(echo $FASTQ | sed 's/.fastq.gz//g')
        bash scripts/01_download.sh \
            $(dirname {output}) \
            $SAMPLE \
            SINGLE \
            {config[GEN_REF]} \
            {config[SRA_CACHE]} \
                &> {log}
        """

rule download_PE:
    output:
        expand(tmp + 'fastq/{sample}/{sample}.fastq_1.gz', zip,
            study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE)
    # params:
        # layout = lambda wildcards:
            # get_metadata(wildcards.sample, layout_col),
    log:
        # base + 'logs/01_download.{sample}.log'
        expand(base + 'logs/01_download.{sample}.log', zip,
            study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE)
    shell:
        """
        FASTQ=$(basename {output})
        SAMPLE=$(echo $FASTQ | sed 's/.fastq_1.gz//g')
        bash scripts/01_download.sh \
            $(dirname {output}) \
            $SAMPLE \
            PAIRED \
            {config[GEN_REF]} \
            {config[SRA_CACHE]} \
                &> {log}
        """

# Rule: expression estimation
rule expression:
    input:
        rules.download_SE.output,
        rules.download_PE.output
    output:
        tmp + 'expression/{sample}/{sample}.quant.sf'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    log:
        base + 'logs/02_expression.{sample}.log'
    shell:
        """
        bash scripts/02_expression.sh \
                 $(dirname {input}) \
                 $(dirname {output}) \
                 {wildcards.sample} \
                 {params.layout} \
                 {config[SALMON_REF]} \
                    &> {log}
        """

# Rule: first-pass alignment
rule align_pass1:
    input:
        rules.download_SE.output,
        rules.download_PE.output
    output:
        tmp + 'junctions/{sample}.junctions.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col)
    log:
        base + 'logs/03_alignment_pass1.{sample}.log'
    shell:
        """
        bash scripts/03_alignment_pass1.sh \
            $(dirname {input}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {params.layout} \
            {config[STAR_REF]} \
               &> {log}
        """

# Rule: second-pass alignment
rule align_pass2:
    input:
        fastq = [rules.download_SE.output, rules.download_PE.output],
        junctions = rules.align_pass1.output
    output:
        tmp + 'alignment/{sample}.bam.tmp'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col),
    log:
        base + 'logs/04_alignment_pass2.{sample}.log'
    shell:
        """
        bash scripts/04_alignment_pass2.sh \
            $(dirname {input.fastq}) \
            $(dirname {input.junctions}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {params.layout} \
            {params.group} \
            {config[STAR_REF]} \
                &> {log}
        """

# Rule: fastq cleanup
rule fastq_cleanup:
    input:
        expression = tmp + 'expression/{sample}/{sample}.quant.sf',
        alignment = tmp + 'alignment/{sample}.bam.tmp'
    output:
        expression = base + 'expression/{sample}/{sample}.quant.sf',
        alignment = tmp + 'alignment/{sample}.bam'
    shell:
        """
        mv $(dirname {input.expression}) $(dirname {output.expression})
        mv {input.alignment} {output.alignment}
        FASTQDIR="{tmp}/data/{wildcards.study}/{wildcards.group}/"
        FASTQDIR+="fastq/{wildcards.sample}"
        rm -r $FASTQDIR
        """

# Rule: variant calling
rule variant_calling:
    input:
        rules.fastq_cleanup.output.alignment
    output:
        tmp + 'variants/{sample}.vcf'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col),
    log:
        base + 'logs/05_variant_calling.{sample}.log'
    shell:
        """
        bash scripts/05_variant_calling.sh \
            $(dirname {input}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {wildcards.group} \
            {config[GEN_REF]} \
            {config[PICARD]} \
            {config[GATK]} \
            {config[KNOWNSNPS]} \
            {config[KNOWNINDELS]} \
                &> {log}
        """

# Rule: alignment cleanup
rule alignment_cleanup:
    input:
        tmp + 'variants/{sample}.vcf'
    output:
        base + 'variants/{sample}.vcf'
    shell:
        """
        mv {input} {output}
        BAMFILE="{tmp}/data/{wildcards.study}/{wildcards.group}/"
        BAMFILE+="alignment/{wildcards.sample}.bam"
        rm $BAMFILE
        """
