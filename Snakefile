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
LAYOUT = config["LAYOUT"]
STUDIES = metadata.loc[metadata[layout_col] == LAYOUT][study_col]
GROUPS = metadata.loc[metadata[layout_col] == LAYOUT][group_col]
SAMPLES = metadata.loc[metadata[layout_col] == LAYOUT][sample_col]

# Set current layout and resulting fastq-type
if LAYOUT == "SINGLE":
    fastq_type = ".fastq.gz"
else:
    fastq_type = "_1.fastq.gz"

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Define directories for finished files / files under processing
finished = 'data/{study}/{group}/'
processingdir = '.processing_snakemake'
processing = processingdir + '/' + finished

# Rule: final outputs
if config["perform_variant_calling"]:
    rule all:
        input:
            expand(finished + 'expression/{sample}/{sample}.quant.sf', zip,
                study = STUDIES, group = GROUPS, sample = SAMPLES),
            expand(finished + 'variants/{sample}.vcf', zip,
                study = STUDIES, group = GROUPS, sample = SAMPLES)
        shell:
            'rm -rf {processingdir}'
else:
    rule all:
        input:
            expand(finished + 'expression/{sample}/{sample}.quant.sf', zip,
                study = STUDIES, group = GROUPS, sample = SAMPLES)
        shell:
            'rm -rf {processingdir}'
    

# Rule: clean
rule clean:
    shell:
        'rm -rf data {processingdir}'

# Rule: download raw data
rule download:
    output:
        processing + 'fastq/{sample}/{sample}' + fastq_type
    log:
        finished + 'logs/01_download.{sample}.log'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    shell:
        """
        bash scripts/01_download.sh \
            $(dirname {output}) \
            {wildcards.sample} \
            {params.layout} \
            {config[GEN_REF]} \
            {config[SRA_CACHE]} \
                %> {log}
        """

# Rule: expression estimation
rule expression:
    input:
        rules.download.output
    output:
        processing + 'expression/{sample}/{sample}.quant.sf'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    log:
        finished + 'logs/02_expression.{sample}.log'
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
rule alignment_pass1:
    input:
        rules.download.output
    output:
        processing + 'junctions/{sample}.junctions.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col)
    log:
        finished + 'logs/03_alignment_pass1.{sample}.log'
    shell:
        """
        bash scripts/03_alignment_pass1.sh \
            $(dirname {input}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {params.layout} \
            {config[STAR_REF]} \
            {config[STAR_THREADS]} \
               &> {log}
        """

# Rule: second-pass alignment
rule alignment_pass2:
    input:
        fastq = rules.download.output,
        junctions = rules.alignment_pass1.output
    output:
        processing + 'alignment/{sample}.bam.tmp'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col),
    log:
        finished + 'logs/04_alignment_pass2.{sample}.log'
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
            {config[STAR_THREADS]} \
                &> {log}
        """

# Rule: fastq cleanup and finalise expression output
if config["perform_variant_calling"]:
    rule finalise_alignment:
        input:
            fastq = rules.download.output,
            expression = rules.expression.output,
            alignment = rules.alignment_pass2.output
        output:
            expression = finished + 'expression/{sample}/{sample}.quant.sf',
            alignment = processing + 'alignment/{sample}.bam'
        shell:
            """
            mv $(dirname {input.expression})/* $(dirname {output.expression})
            mv {input.alignment} {output.alignment}
            rm -r $(dirname {input.fastq})
            """
else:
    rule finalise_expression:
        input:
            fastq = rules.download.output,
            expression = rules.expression.output
        output:
            finished + 'expression/{sample}/{sample}.quant.sf',
        shell:
            """
            mv $(dirname {input.expression})/* $(dirname {output})
            rm -r $(dirname {input.fastq})
            """
    
# Rule: variant calling
if config["perform_variant_calling"]:
    rule variant_calling:
        input:
            rules.finalise_alignment.output.alignment
        output:
            processing + 'variants/{sample}.vcf'
        params:
            layout = lambda wildcards:
                get_metadata(wildcards.sample, layout_col),
            group = lambda wildcards:
                get_metadata(wildcards.sample, group_col),
        log:
            finished + 'logs/05_variant_calling.{sample}.log'
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

# Rule: alignment cleanup and finalise variant calling output
if config["perform_variant_calling"]:
    rule finalise_variant_calling:
        input:
            alignment = rules.finalise_alignment.output.alignment,
            variants = rules.variant_calling.output
        output:
            finished + 'variants/{sample}.vcf'
        shell:
            """
            mv {input.variants} {output}
            rm -r $(dirname {input.alignment})
            """
