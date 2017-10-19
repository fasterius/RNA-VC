#!/usr/bin/env python3

# Packages
import pandas as pd

# Mode
mode = 'SRR'

# Metadata
metadata_file = 'test.metadata.txt'
metadata = pd.read_table(metadata_file, encoding = 'ISO-8859-1')

# Get SRR and GSEs
SRR = metadata['SRR']
SRR_SE = metadata.loc[metadata['LibraryLayout'] == 'SINGLE']['SRR']
SRR_PE = metadata.loc[metadata['LibraryLayout'] == 'PAIRED']['SRR']
GSE_SE = metadata.loc[metadata['LibraryLayout'] == 'SINGLE']['GSE']
GSE_PE = metadata.loc[metadata['LibraryLayout'] == 'PAIRED']['GSE']

rule all:
    input:
        expand('analysis/variants/{SRR}/{SRR}.vcf', SRR = SRR)

# Download and estimate expression
rule download:
    output:
        expand('data/fastq_se/{SRR}/{SRR}.fastq.gz',
               SRR = SRR_SE),
        expand('data/fastq_pe/{SRR}/{SRR}_1.fastq.gz',
               SRR = SRR_PE),
        expand('analysis/expression/{SRR}/{SRR}.abundance.tsv',
               SRR = SRR)
    shell:
        "touch {output}"
        # "scripts/01_download_fastq.sh {wildcards.sample} SINGLE"

# Alignment
rule align:
    input:
        expand('data/fastq_se/{SRR}/{SRR}.fastq.gz',
               SRR = SRR_SE),
        expand('data/fastq_pe/{SRR}/{SRR}_1.fastq.gz',
               SRR = SRR_PE)
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
