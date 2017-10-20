#!/usr/bin/env python3

# Packages
import pandas as pd

# Config file
configfile: "config.yaml"

# Read metadata
metadata_file = config["metadata"]
metadata = pd.read_table(metadata_file, encoding = 'ISO-8859-1')

# Get relevant column names from config
sample_col = config["sample_col"]
group_col = config["group_col"]
layout_col = config["layout_col"]

# Get samples and groups from metadata
SAMPLES = metadata[sample_col]
SAMPLES_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][sample_col]
SAMPLES_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][sample_col]
GROUPS_SE = metadata.loc[metadata[layout_col] == 'SINGLE'][group_col]
GROUPS_PE = metadata.loc[metadata[layout_col] == 'PAIRED'][group_col]

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Rules

# Rule: create all outputs
rule all:
    input:
        expand('data/fastq/{sample_SE}/{sample_SE}.fastq.gz',
            sample_SE = SAMPLES_SE),
        expand('data/fastq/{sample_PE}/{sample_PE}_1.fastq.gz',
            sample_PE = SAMPLES_PE),
        expand('results/{sample}.vcf', sample = SAMPLES)
    
# Rule: download and estimate expression on single-end data
rule download_SE:
    output:
        'data/fastq/{sample_SE}/{sample_SE}.fastq.gz',
        'data/expression/{sample_SE}/{sample_SE}.abundance.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample_SE, layout_col)
    shell:
        "touch {output}"
        # """
        # bash scripts/01_download_fastq.sh {wildcards.sample_SE} \
            # {params.layout} {config[GEN_REF]} {config[SRA_CACHE]}
        # """

# Rule: download and estimate expression on paired-end data
rule download_PE:
    output:
        'data/fastq/{sample_PE}/{sample_PE}_1.fastq.gz',
        'data/expression/{sample_PE}/{sample_PE}.abundance.tsv'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample_PE, layout_col)
    shell:
        "touch {output}"
        # """
        # bash scripts/01_download_fastq.sh {wildcards.sample_PE} \
            # {params.layout} {config[GEN_REF]} {config[SRA_CACHE]}
        # """

# Rule: single-end alignment
rule align_SE:
    input:
        'data/fastq/{sample_SE}/{sample_SE}(?!_1).fastq.gz'
    output:
        'data/alignment/{sample_SE}/{sample_SE}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample_SE, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample_SE, group_col)
    shell:
        'touch {output}'
        # """
        # bash scripts/02_alignment.sh {wildcards.sample_SE} \
            # {params.layout} {params.group} \
            # {config[GEN_REF]} {config[STAR_REF]}
        # """

# Rule: paired-end alignment
rule align_PE:
    input:
        'data/fastq/{sample_PE}/{sample_PE}_1.fastq.gz'
    output:
        'data/alignment/{sample_PE}/{sample_PE}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample_PE, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample_PE, group_col)
    shell:
        'touch {output}'
        # """
        # bash scripts/02_alignment.sh {wildcards.sample_PE} \
            # {params.layout} {params.group} \
            # {config[GEN_REF]} {config[STAR_REF]}
        # """

# Rule: variant calling
rule variant_calling:
    input:
        'data/alignment/{sample}/{sample}.bam'
    output:
        'results/{sample}.vcf'
    shell:
        'touch {output}'
