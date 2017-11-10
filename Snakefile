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
sample_col = config["sample_col"]
layout_col = config["layout_col"]

# Get studies and samples from metadata
LAYOUT = config["LAYOUT"]
STUDIES = metadata.loc[metadata[layout_col] == LAYOUT][study_col]
SAMPLES = metadata.loc[metadata[layout_col] == LAYOUT][sample_col]

# Function for fetching specified metadata from a sample
def get_metadata(sample, column):
    result = metadata.loc[metadata[sample_col] == sample][column].values[0]
    return result

# Set current layout and resulting fastq-type
if LAYOUT == "SINGLE":
    fastq_type = ".fastq.gz"
else:
    fastq_type = "_1.fastq.gz"

# Define directories for outdir files / files under outdir
outdir = 'data/{study}/'
tempdir = '.snakemake_temp/' + outdir

# Rule: final outputs
if config["perform_variant_calling"]:
    rule all:
        input:
            expand(outdir + 'expression/{sample}/{sample}.quant.sf', zip,
                study = STUDIES, sample = SAMPLES),
            expand(outdir + 'variants/{sample}.vcf.gz', zip,
                study = STUDIES, sample = SAMPLES)
        shell:
            "rm -r {tempdir}"
else:
    rule all:
        input:
            expand(outdir + 'expression/{sample}/{sample}.quant.sf', zip,
                study = STUDIES, sample = SAMPLES)
        shell:
            "rm -r {tempdir}"

# Rule: clean
rule clean:
    shell:
        'rm -rf data {tempdir}'

# Rule: download raw data
rule download:
    output:
        fastq = tempdir + 'fastq/{sample}/{sample}' + fastq_type,
        expression = outdir + 'expression/{sample}/{sample}.quant.sf'
    log:
        outdir + 'logs/{sample}.01_download.log'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    shell:
        """
        bash scripts/01_download.sh \
            $(dirname {output.fastq}) \
            $(dirname {output.expression}) \
            {wildcards.sample} \
            {params.layout} \
            {config[SALMON_REF]} \
            {config[SRA_CACHE]} \
                &> {log}
        """

# Rule: alignment
rule alignment:
    input:
        rules.download.output.fastq
    output:
        tempdir + 'alignment/{sample}/{sample}.bam'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    log:
        outdir + 'logs/{sample}.02_alignment.log'
    shell:
        """
        bash scripts/02_alignment.sh \
            $(dirname {input}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {params.layout} \
            {config[STAR_REF]} \
            {config[STAR_THREADS]} \
            {log} \
               &> {log}
        """

# Rule: variant calling
rule variant_calling:
    input:
        rules.alignment.output
    output:
        outdir + 'variants/{sample}.vcf.gz'
    params:
        layout = lambda wildcards:
            get_metadata(wildcards.sample, layout_col),
    log:
        outdir + 'logs/{sample}.03_variant_calling.log'
    shell:
        """
        bash scripts/03_variant_calling.sh \
            $(dirname {input}) \
            $(dirname {output}) \
            {wildcards.sample} \
            {config[GEN_REF]} \
            {config[PICARD]} \
            {config[GATK]} \
            {config[KNOWNSNPS]} \
            {config[KNOWNINDELS]} \
            {config[SNPEFF]} \
            {config[SNPEFFASSEMBLY]} \
                &> {log}
        """
