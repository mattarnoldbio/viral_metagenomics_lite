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

**mattarnoldbio/viral_metagenomics_lite** is a lightweight viral metagenomics workflow for short-read sequencing data. It performs read QC and trimming, removes host reads, assembles contigs, searches for viral matches against reference databases, assigns taxonomy, and filters the final hit tables into a clean set of outputs for downstream interpretation.

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
- A DIAMOND database file for protein search (`--diamond_db`)
- A Krona database for taxonomy plots (`--krona_db`)
- A BLAST nucleotide database basename/prefix for confirmatory nt search (`--blast_db`)
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
- `ref_genome`: Bowtie2 index prefix used for host depletion, not the full path to a single `.bt2` file

The pipeline will auto-detect paired-end versus single-end reads from the file information in the samplesheet.

For host depletion, give the common prefix used to build the Bowtie2 index, e.g. `/path/to/host_index/human_hg38`. The pipeline looks for files matching `${ref_genome}*.bt2`, so values like `/path/to/host_index/human_hg38` are correct, while a full filename such as `/path/to/host_index/human_hg38.1.bt2` would not be appropriate.

The BLAST database is also specified as a basename/prefix rather than a single file. For example, if your nt database exists as `/path/to/nt.nhr`, `/path/to/nt.nin`, `/path/to/nt.nsq`, etc., then pass `--blast_db /path/to/nt` rather than a full file path like `/path/to/nt.nsq`. In contrast, `--diamond_db` should point to the actual database file, usually ending in `.dmnd`.

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
- `--diamond_db`: full path to the DIAMOND database file, usually ending in `.dmnd`
- `--krona_db`: path to the Krona taxonomy database used to assign taxonomic labels
- `--blast_db`: BLAST database prefix/basename, e.g. `/path/to/nt` for an nt database whose files are `/path/to/nt.nhr`, `/path/to/nt.nin`, `/path/to/nt.nsq`, etc.
- `--min_contig_length`: minimum contig length retained after filtering

Example with a params file:

```bash
nextflow run mattarnoldbio/viral_metagenomics_lite \
  -profile docker \
  -params-file params.yaml
```

A complete example parameter file is provided in [params.yaml](params.yaml). It contains the default values from the pipeline config and can be copied and edited for your own run.

```yaml
input: './samplesheet.csv'
outdir: './results'
diamond_db: '/path/to/nr.dmnd'
krona_db: '/path/to/krona_db'
blast_db: '/path/to/nt'
min_contig_length: 1000
```

> [!IMPORTANT]
> In newer Nextflow versions, non-string CLI values are not always parsed reliably from the command line. In practice, settings such as `--local_max_cpus`, `--local_max_memory`, `--min_contig_length`, and other numeric or typed values are usually best provided in a YAML/JSON params file rather than passed directly on the CLI.
>
> For example, if you need to set the minimum contig length or memory/thread limits, prefer:
>
> ```bash
> nextflow run mattarnoldbio/viral_metagenomics_lite \
>   -profile docker \
>   -params-file params.yaml
> ```
>
> with entries such as:
>
> ```yaml
> local_max_cpus: 10
> local_max_memory: '24.GB'
> min_contig_length: 1000
> ```
>
> This avoids problems caused by the CLI treating some values as strings or otherwise failing to coerce them into the expected types.

### Output publishing modes

The pipeline exposes Nextflow's `publish_dir_mode` option, which controls how output files are written to the output directory. This can be set on the command line as:

```bash
nextflow run mattarnoldbio/viral_metagenomics_lite \
  -profile docker \
  --input ./samplesheet.csv \
  --outdir ./results \
  --publish_dir_mode symlink
```

Common choices are:

- `copy`: default behavior; duplicates files into the output directory. Safest and most portable when you want fully independent output files.
- `symlink`: creates symbolic links to files instead of copying them. Useful when disk space is limited and the outputs are still being used from the original working directory.
- `rellink`: creates relative symlinks to the original files. Helpful when you want a lightweight link while keeping the output directory portable across a project tree.
- `link`: creates hard links or symlinks depending on the filesystem and setup; less common and best used when you know the backend supports it.
- `move`: moves files into the output directory rather than copying them. Appropriate when you want to clean up intermediates or keep only the published result set.
- `copyNoFollow`: copies files but does not follow symlinks. Useful when you want a standard copy without dereferencing link targets.

In practice, `copy` is the safest default for most users, while `symlink` or `rellink` are convenient when storage space matters and the files are not expected to be moved elsewhere. `move` is more aggressive and should be used with care because it changes file location rather than leaving the original in the work directory.

### `-resume` and publish mode

Nextflow's `-resume` flag reuses cached results when a previous run has the same inputs and parameters. This is especially useful when re-running a pipeline after a small parameter change or after fixing a failed step.

The important point for this pipeline is that `-resume` helps reuse computation, but it does not automatically change the way files are published. If you change `--publish_dir_mode` between runs, you may still get the cached downstream results, while the output directory contents are republished according to the new mode. In particular:

- `copy` will write independent copies into the output directory
- `symlink`/`rellink` will publish lightweight links, which can be convenient but may look stale or disappear if the original work files are cleaned up
- `move` can make the published output feel more "final" but is the most disruptive choice when resuming or rerunning a workflow

A common pattern is to use `-resume` for iterative debugging or reruns, while keeping `--publish_dir_mode copy` for reproducible, portable outputs.

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
