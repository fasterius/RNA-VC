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

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Wildcard constraints
wildcard_constraints:
    # fastq_SE = ".+_0",
    # fastq_PE = ".+_1",
    # sample_SE = ".+(?!_1)",
    # sample_PE = ".+(?!_0)"

# Rule: create all outputs
rule all:
    input:
        # expand('data/{study}/fastq/{sample}/{sample}_0.fastq.gz', zip,
            # study = STUDIES_SE, sample = SAMPLES_SE),
        # expand('data/{study}/fastq/{sample}/{sample}_1.fastq.gz', zip,
            # study = STUDIES_PE, sample = SAMPLES_PE),
        expand('data/se/{study}/variants/{sample}.vcf', zip,
            study = STUDIES_SE, sample = SAMPLES_SE),
        expand('data/pe/{study}/variants/{sample}.vcf', zip,
            study = STUDIES_PE, sample = SAMPLES_PE)
    
# Rule: download and estimate expression on single-end data
rule download_SE:
    output:
        '.tmp/se/data/{study}/fastq/{sample}/{sample}_0.fastq.gz'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col)
    shell:
        "touch {output}"
        # """
        # bash scripts/01_download_fastq.sh {wildcards.sample} \
            # {params.layout} {config[GEN_REF]} {config[SRA_CACHE]}
        # """

rule download_PE:
    output:
        '.tmp/pe/data/{study}/fastq/{sample}/{sample}_1.fastq.gz'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col)
    shell:
        "touch {output}"
        # """
        # bash scripts/01_download_fastq.sh {wildcards.sample} \
            # {params.layout} {config[GEN_REF]} {config[SRA_CACHE]}
        # """

# Rule: single-end alignment (pass 2)
rule align_SE_pass1:
    input:
        rules.download_SE.output
    output:
        '.tmp/se/data/{study}/alignment/{sample}/{sample}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col)
    shell:
        'touch {output}'
        # """
        # bash scripts/02_alignment.sh {wildcards.sample} \
            # {params.layout} {params.group} \
            # {config[GEN_REF]} {config[STAR_REF]}
        # """

# Rule: paired-end alignment (pass 1)
rule align_PE_pass1:
    input:
        rules.download_PE.output
    output:
        '.tmp/pe/data/{study}/alignment/{sample}/{sample}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
        group = lambda wildcards:
            get_metadata(wildcards.sample, group_col)
    shell:
        'touch {output}'
        # """
        # bash scripts/02_alignment.sh {wildcards.sample_PE} \
            # {params.layout} {params.group} \
            # {config[GEN_REF]} {config[STAR_REF]}
        # """

# Rule: variant calling
rule variant_calling_SE:
    input:
        rules.align_SE_pass1.output
    output:
        '.tmp/se/data/{study}/variants/{sample}.vcf'
    shell:
        'touch {output}'

# Rule: variant calling
rule variant_calling_PE:
    input:
        rules.align_PE_pass1.output
    output:
        '.tmp/pe/data/{study}/variants/{sample}.vcf'
    shell:
        'touch {output}'

rule rename_SE:
    input:
        '.tmp/se/data/{study}'
        # '.tmp/se/data/{study}/fastq/{sample_SE}/{sample_SE}_0.fastq.gz'
    output:
        'data/se/{study}'
        # 'data/{study}/fastq/{sample_SE}/{sample_SE}_0.fastq.gz'
    shell:
        "mv {input} {output}"

rule rename_PE:
    input:
        '.tmp/pe/data/{study}'
        # '.tmp/pe/data/{study}/fastq/{sample_PE}/{sample_PE}_1.fastq.gz'
    output:
        'data/pe/{study}'
        # 'data/{study}/fastq/{sample_PE}/{sample_PE}_1.fastq.gz'
    shell:
        "mv {input} {output}"
