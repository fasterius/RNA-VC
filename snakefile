#!/usr/bin/env python3

# Packages
import pandas as pd

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

# Get samples and groups from metadata
STUDIES = metadata[study_col]
STUDIES_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][study_col]
STUDIES_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][study_col]
GROUPS = metadata[group_col]
GROUPS_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][group_col]
GROUPS_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][group_col]
SAMPLES = metadata[sample_col]
SAMPLES_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][sample_col]
SAMPLES_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][sample_col]

# # Function for specifying layout-dependent file names
# def input_sample_layouts(wildcards):
#
    # # Get current sample
    # sample = wildcards.sample
    # study = wildcards.study
#
    # # Get sample layout
    # layout = metadata.loc[metadata[sample_col] == sample][layout_col].values[0]
#
    # # # Get sample parameters
    # # study = metadata.loc[metadata[sample_col] == sample][study_col].values[0]
#
    # # Sample string
    # string = '.tmp/data/{study}/fastq/{sample}/{sample}_1.fastq.gz'
    # string = string.format(sample=sample, study=study)
#
    # # Get sample names
    # if layout == 'SINGLE':
        # samples = [string.replace('_1', '')]
    # elif layout == 'PAIRED':
        # samples = [string, string.replace('_1', '_2')]
    #
    # # Return sample names
    # return samples

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Strings for base directories
base = 'data/{study}/{group}/'
tmp = '.tmp/' + base

# Rule: final outputs
rule all:
    input:
        expand(base + 'expression/{sample}.abundance.tsv', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES),
        expand(base + 'variants/{sample}.vcf', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES)

# Rule: clean
rule clean:
    shell:
        'rm -r data .tmp'

# Rule: download and estimate expression on single-end data
rule download:
    output:
        expand(tmp + 'fastq/{sample}.fastq.gz', zip,
            study = STUDIES_SE, group = GROUPS_SE, sample = SAMPLES_SE),
        expand(tmp + 'fastq/{sample}_1.fastq.gz', zip,
            study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE),
        expand(tmp + 'fastq/{sample}_2.fastq.gz', zip,
            study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE),
        expand(base + 'expression/{sample}.abundance.tsv', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES)
    log:
        expand(base + 'logs/download.{sample}.log', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES)
    run:
        for current_output in output:
            sample = current_output.split('/')[3]
            shell('touch {current_output} 2> {log}')
        # """
        # bash scripts/01_download_fastq.sh {wildcards.sample} \
            # {params.layout} {config[GEN_REF]} {config[SRA_CACHE]}
        # """

# Rule: first-pass alignment
rule align_pass1:
    input:
        rules.download.output
    output:
        tmp + 'junctions/{sample}.junctions.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col)
    log:
        base + 'logs/alignment_pass1.{sample}.log'
    shell:
        'touch {output} 2> {log}'
        # """
        # bash scripts/02_alignment.sh {wildcards.sample} \
            # {params.layout} {params.group} \
            # {config[GEN_REF]} {config[STAR_REF]}
        # """

# Rule: second-pass alignment
rule align_pass2:
    input:
        rules.download.output,
        rules.align_pass1.output
    output:
        tmp + 'alignment/{sample}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col),
    log:
        base + 'logs/alignment_pass2.{sample}.log'
    shell:
        'touch {output} 2> {log}'

# Rule: variant calling
rule variant_calling:
    input:
        rules.align_pass2.output
    output:
        base + 'variants/{sample}.vcf'
    log:
        base + 'logs/variant_calling.{sample}.log'
    shell:
        'touch {output} 2> {log}'
