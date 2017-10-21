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
# def input_sample_layout(wildcards):
#
    # # Get current sample
    # sample = wildcards.sample
    # study = wildcards.study
#
    # # Get sample layout
    # layout = metadata.loc[metadata[sample_col] == sample][layout_col].values[0]
#
    # # Sample string
    # string = '.tmp/data/{study}/fastq/{sample}_1.fastq.gz'
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
            study = STUDIES_PE, group = GROUPS_PE, sample = SAMPLES_PE)
    log:
        expand(base + 'logs/download.{sample}.log', zip,
            study = STUDIES, group = GROUPS, sample = SAMPLES)
    run:
        for path in output:

            # Define file path
            separated_path = path.split('/')

            # Get current sample
            fastq = separated_path[-1]
            sample = re.sub('_[1-2]', '', fastq.split('.')[0])

            # Get file path
            path = path.replace('/' + fastq, '')

            # Get current layout
            current = metadata.loc[metadata[sample_col] == sample]
            layout = current[layout_col].values[0]
            
            # Execute download script with current parameters
            shell('touch {output}')
            # shell("bash scripts/01_download.sh \
                    # {path} \
                    # {sample} \
                    # {layout} \
                    # {config[GEN_REF]} \
                    # {config[SRA_CACHE]} \
                        # 2>&1 {log}")

# Rule: expression estimation
rule expression:
    input:
        rules.download.output
    output:
        base + 'expression/{sample}.abundance.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    log:
        base + 'logs/expression.{sample}.log'
    shell:
        'touch {output}'
        # """
        # bash scripts/02_expression.sh \
                 # $(dirname {output}) \
                 # {wildcards.sample} \
                 # {params.layout} \
                 # {config[GEN_REF]} \
                    # 2>&1 {log}
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
        'touch {output}'
        # """
        # bash scripts/03_alignment_pass1.sh \
            # $(dirname {output}) \
            # {wildcards.sample} \
            # {params.layout} \
            # {params.group} \
            # {config[GEN_REF]} \
            # {config[STAR_REF]} \
               # 2>&1 {log}
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
