# RNA-VC
[![License: MIT][badge]][licence]

## Overview

**RNA-VC** is a pipeline that can analyse any RNA sequencing (RNA-seq) raw data
available in the *Gene Expression Omnibus* ([GEO][geo]), yielding both
variant calling data and gene/transcript expression data. It can also be run
without the variant calling, if you are only interested in expression data.
It is written with a combination of `bash` and `python`, all wrapped up in the
`snakemake` workflow manager system.

RNA-VC differentiates itself from the other available RNA-seq variant calling
pipelines in that you do not need to have the raw data downloaded before you
start; RNA-VC takes care of that for you. All you need is a metadata file with
the specified samples you want to analyse.

RNA-VC is a *work in progress*, and does not come with any guarantees. I am
using this pipeline in my own work, but just because it runs successfully on
my system does not mean it will do so on yours. Nevertheless, I share it here
in the hope that it may prove useful for others.

## Usage

There are a number of bioinformatic software packages that need to be installed
in order to run the pipeline:

 * [Python 3][python] and [Snakemake][snakemake] for running the pipeline
 * [SRAtools][sratools] for downloading raw fastq files
 * [Salmon][salmon] for estimating expression level
 * [STAR][star] for aligning the reads
 * [SAMtools][samtools] for working with the resulting alignment files
 * [PICARD][picard] for marking read duplicates
 * [GATK][gatk] for performing variant calling
 * [snpEFF][snpeff] for annotating resulting variants

You must have all of these installed if you are to run the pipeline (if you are
using a computer cluster they might already be installed, if you're lucky).
After you have made sure they are all installed, you can `clone` this
repository to get all the relevant files:

```bash
git clone https://github.com/fasterius/RNA-VC
```

You also need to provide RNA-VC with the metadata describing the GEO samples
you want to analyse. This means that you need to provide, at the very least,
*SRR IDs*, their corresponding *GSE IDs*, their read layouts (listing them 
either as "SINGLE" or "PAIRED", respectively) and, lastly, any subgroup the
samples belong to. Such a metadata file might look like this:

| Study     | Sample     | Group     | Layout     |
| --------- | ---------- | --------- | ---------- |
| GSE81194  | SRR3479755 | HCT116    | PAIRED     |
| GSE81194  | SRR3479758 | RKO       | PAIRED     |
| GSE81194  | SRR3479759 | RKO       | PAIRED     |

The `config.yaml` file provided can then be edited according to the structure
of your metadata file, in addition to the locations of the references, indexes
software paths required. You finally run the pipeline using `snakemake`:

```bash
snakemake --config LAYOUT=PAIRED
```

If you have mixed read layouts you have to run the pipeline twice: once for
each read layout. Alternatively, if you are running the pipeline on a cluster
(RNA-VC currently only supports `SLURM`) you can use the `submit_snakemake.sh`
wrapper, which will automatically loop through both read layouts:

```bash
bash submit_snakemake.sh
```

## License

The pipeline is available with a MIT licence. It is free software: you may
redistribute it and/or modify it under the terms of the MIT license. For more
information, please see the `LICENCE` file.

[badge]: https://img.shields.io/badge/license-mit-blue.svg
[licence]: https://opensource.org/licenses/mit
[geo]: https://www.ncbi.nlm.nih.gov/geo/
[python]: https://www.python.org/
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[sratools]: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
[salmon]: https://combine-lab.github.io/salmon/
[star]: https://github.com/alexdobin/STAR
[samtools]: http://www.htslib.org/
[picard]: http://broadinstitute.github.io/picard/
[gatk]: https://software.broadinstitute.org/gatk/
[snpeff]: http://snpeff.sourceforge.net/
