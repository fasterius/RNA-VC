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

rule all:
    input:
        expand('data/{SRR_SE}/{SRR_SE}.fastq.gz', SRR_SE = SRR_SE),
        expand('data/{SRR_PE}/{SRR_PE}_1.fastq.gz', SRR_PE = SRR_PE),
        expand('analysis/variants/{SRR}/{SRR}.vcf', SRR = SRR)

# Function for matching metadata to current sample
def get_metadata(sample, column):
    result = metadata.loc[metadata["SRR"] == sample][column].values
    return result

rule download_PE:
    output:
        'data/{SRR_PE}/{SRR_PE}_1.fastq.gz',
        'analysis/expression/{SRR_PE}/{SRR_PE}.abundance.tsv'
    params:
        layout = lambda wildcards: get_metadata(wildcards.SRR_PE,
                                                "LibraryLayout")
    shell:
        "bash testscript.sh {wildcards.SRR_PE} {params.layout}"
    
# Download and estimate expression
rule download_SE:
    output:
        'data/{SRR_SE}/{SRR_SE}.fastq.gz',
        'analysis/expression/{SRR_SE}/{SRR_SE}.abundance.tsv'
    params:
        layout = lambda wildcards: get_metadata(wildcards.SRR_SE,
                                                "LibraryLayout")
    shell:
        "bash testscript.sh {wildcards.SRR_SE} {params.layout}"

# Alignment
rule align:
    input:
        expand('data/{SRR}/{SRR}.fastq.gz', SRR = SRR_SE),
        expand('data/{SRR}/{SRR}_1.fastq.gz', SRR = SRR_PE)
    output:
        expand('analysis/alignment/{SRR}/{SRR}.bam', SRR = SRR)
    shell:
        'touch {output}'

# Variant calling
rule variant_calling:
    input:
        expand('analysis/alignment/{SRR}/{SRR}.bam', SRR = SRR)
    output:
        expand('analysis/variants/{SRR}/{SRR}.vcf', SRR = SRR)
    shell:
        'touch {output}'
