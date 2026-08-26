# mattarnoldbio/viral_metagenomics_lite

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/mattarnoldbio/viral_metagenomics_lite)
[![GitHub Actions CI Status](https://github.com/mattarnoldbio/viral_metagenomics_lite/actions/workflows/nf-test.yml/badge.svg)](https://github.com/mattarnoldbio/viral_metagenomics_lite/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/mattarnoldbio/viral_metagenomics_lite/actions/workflows/linting.yml/badge.svg)](https://github.com/mattarnoldbio/viral_metagenomics_lite/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/mattarnoldbio/viral_metagenomics_lite)

## Introduction

**mattarnoldbio/viral_metagenomics_lite** is a lightweight metagenomics workflow for identifying viral contigs from short-read sequencing data. It takes FASTQ samples, runs read quality checks and trimming, removes host reads using a Bowtie2 index, assembles contigs with MEGAHIT, searches contigs against protein and nucleotide databases, classifies hits with Krona, and filters the final set of viral and non-viral hit calls into a clean output directory.

The pipeline is designed for a typical viromics or virus-discovery workflow in which the user provides a samplesheet of raw read files together with host depletion references and pre-built sequence databases. It produces QC summaries, trimmed/deduplicated reads, host-depleted reads, assembled contigs, taxonomy plots, and filtered hit tables for downstream interpretation.

### Default workflow steps

1. Read QC with FastQC
2. Adapter trimming and quality filtering with Trim Galore
3. Duplicate/artifact filtering with PRINSEQ++
4. Host read depletion using Bowtie2 against a per-sample reference index
5. Assembly of cleaned reads into contigs with MEGAHIT
6. Protein homology search against a DIAMOND database
7. Nucleotide search against a BLAST database
8. Taxonomic assignment and Krona plot generation
9. Virus-hit extraction and contig filtering
10. MultiQC summary report generation

## Requirements

Before running the pipeline, you need:

- A samplesheet listing sample IDs and FASTQ file locations
- A Bowtie2 host-reference index prefix for host depletion (`ref_genome` in the samplesheet)
- A DIAMOND database for protein search (`--diamond_db`)
- A Krona database for taxonomy plots (`--krona_db`)
- A BLAST nucleotide database for confirmatory nt search (`--blast_db`)
- An output directory (`--outdir`)

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set up Nextflow. Make sure to test your setup before running the workflow on real data.

## Samplesheet format

Create a CSV samplesheet with at least these columns:

```csv
sample,fastq_1,fastq_2,ref_genome
SAMPLE_01,/path/to/sample_01_R1.fastq.gz,/path/to/sample_01_R2.fastq.gz,/path/to/host_reference/index_prefix
SAMPLE_02,/path/to/sample_02_R1.fastq.gz,,/path/to/host_reference/index_prefix
```

Field descriptions:

- `sample`: unique sample identifier
- `fastq_1`: path to read 1 (required)
- `fastq_2`: path to read 2 for paired-end data; leave blank for single-end samples
- `ref_genome`: path prefix for the Bowtie2 index used for host depletion

The pipeline will auto-detect paired-end versus single-end reads from the file information in the samplesheet.

## Running the pipeline

A minimal run looks like this:

```bash
nextflow run mattarnoldbio/viral_metagenomics_lite \
  -profile docker \
  --input ./samplesheet.csv \
  --outdir ./results \
  --diamond_db /path/to/nr.dmnd \
  --krona_db /path/to/krona_db \
  --blast_db /path/to/nt \
  --min_contig_length 1000
```

The most relevant arguments to specify are:

- `--input`: path to the sample sheet CSV
- `--outdir`: results output directory
- `--diamond_db`: DIAMOND database used for BLASTX-style homology search
- `--krona_db`: Krona taxonomy database used to assign taxonomic labels
- `--blast_db`: BLAST nt database used for secondary nucleotide confirmation
- `--min_contig_length`: minimum contig length retained after filtering

Example with a params file:

```bash
nextflow run mattarnoldbio/viral_metagenomics_lite \
  -profile docker \
  -params-file params.yaml
```

```yaml
input: './samplesheet.csv'
outdir: './results'
diamond_db: '/path/to/nr.dmnd'
krona_db: '/path/to/krona_db'
blast_db: '/path/to/nt'
min_contig_length: 1000
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters; see the nf-core documentation for details.

## Output structure

The pipeline writes results into the output directory, including:

- `multiqc/` and a MultiQC HTML report
- `preprocessed_reads/` for trimmed/deduplicated and host-depleted reads
- `contigs/` for separate and combined assembled contigs, including viral contigs
- `db_search_results/` for DIAMOND and BLAST output tables
- `final_results/` for the filtered final virus/non-virus hit sets
- `pipeline_info/` for execution metadata and software versions

## Credits

mattarnoldbio/viral_metagenomics_lite was originally written by Matt Arnold.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used in the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
