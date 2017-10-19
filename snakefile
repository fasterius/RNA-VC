#!/usr/bin/env python3

# Packages
import pandas as pd

# Config file
configfile: "config.yaml"

# Read metadata
metadata_file = config["metadata"]
metadata = pd.read_table(metadata_file, encoding = 'ISO-8859-1')

# Get SRR and GSEs
SRR = metadata['SRR']
SRR_SE = metadata.loc[metadata['LibraryLayout'] == 'SINGLE']['SRR']
SRR_PE = metadata.loc[metadata['LibraryLayout'] == 'PAIRED']['SRR']
GSE_SE = metadata.loc[metadata['LibraryLayout'] == 'SINGLE']['GSE']
GSE_PE = metadata.loc[metadata['LibraryLayout'] == 'PAIRED']['GSE']

# Function for matching metadata to current sample
def get_metadata(sample, column):
    result = metadata.loc[metadata["SRR"] == sample][column].values
    return result

# Rules

# Rule: create all outputs
rule all:
    input:
        expand('data/fastq/{SRR_SE}/{SRR_SE}.fastq.gz', SRR_SE = SRR_SE),
        expand('data/fastq/{SRR_PE}/{SRR_PE}_1.fastq.gz', SRR_PE = SRR_PE),
        expand('results/{SRR}.vcf', SRR = SRR)

# Rule: download and estimate expression on paired-end data
rule download_PE:
    output:
        'data/fastq/{SRR_PE}/{SRR_PE}_1.fastq.gz',
        'data/expression/{SRR_PE}/{SRR_PE}.abundance.tsv'
    params:
        layout = lambda wildcards: get_metadata(wildcards.SRR_PE,
                                                "LibraryLayout")
    shell:
        "bash testscript.sh {wildcards.SRR_PE} {params.layout}"
    
# Rule: download and estimate expression on single-end data
rule download_SE:
    output:
        'data/fastq/{SRR_SE}/{SRR_SE}.fastq.gz',
        'data/expression/{SRR_SE}/{SRR_SE}.abundance.tsv'
    params:
        layout = lambda wildcards: get_metadata(wildcards.SRR_SE,
                                                "LibraryLayout")
    shell:
        "bash testscript.sh {wildcards.SRR_SE} {params.layout}"

# Rule: alignment
rule align:
    input:
        expand('data/fastq/{SRR}/{SRR}.fastq.gz', SRR = SRR_SE),
        expand('data/fastq/{SRR}/{SRR}_1.fastq.gz', SRR = SRR_PE)
    output:
        expand('data/alignment/{SRR}/{SRR}.bam', SRR = SRR)
    shell:
        'touch {output}'

# Rule: variant calling
rule variant_calling:
    input:
        expand('data/alignment/{SRR}/{SRR}.bam', SRR = SRR)
    output:
        expand('results/{SRR}.vcf', SRR = SRR)
    shell:
        'touch {output}'
