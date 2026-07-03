/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE                } from '../modules/nf-core/trimgalore/main'
include { PRINSEQPLUSPLUS           } from '../modules/nf-core/prinseqplusplus/main'
include { BOWTIE2_BUILD             } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN             } from '../modules/nf-core/bowtie2/align/main'
include { MEGAHIT                   } from '../modules/nf-core/megahit/main'
include { CONCATENATE_CONTIG_FILES  } from '../modules/local/concatenate_contig_files/main' 
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_viral_metagenomics_lite_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRAL_METAGENOMICS_LITE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // MODULE: Trim adapters and low quality bases with Trim Galore
    // TODO: Check params
    TRIMGALORE(ch_samplesheet)

    //
    // MODULE: Remove PCR duplicates and low quality reads with PRINSEQ++
    // TODO: Check params
    PRINSEQPLUSPLUS(ch_samplesheet)
    
    ch_samplesheet
        .multiMap { meta, reads ->
            reads_ch: [ meta, reads ]
            index_ch: [
                meta,                              // reuse same meta — fine, since script ignores meta2 content
                file("${meta.ref_genome}*.bt2")     // glob matches all 6 index files for this sample's ref genome
            ]
        }
        .set { ch_for_align }

    // ch_for_align.reads_ch.view()
    // ch_for_align.index_ch.view()

    ch_fasta = channel.of([ [:], [] ]).first() // create a channel with a single tuple of empty meta and empty fasta path, to be used in BOWTIE2_ALIGN

    //
    // MODULE: Exlcude host reads with Bowtie2
    // TODO: Check params
    // TODO: Test with single-end reads
    BOWTIE2_ALIGN(ch_for_align.reads_ch, ch_for_align.index_ch, ch_fasta, true, true)

    ch_reads_host_depleted = BOWTIE2_ALIGN.out.fastq.map {meta, reads ->
        def r1 = null
        def r2 = null

        if (meta.single_end) {
            r1 = reads
        } else {
            r1 = reads[0]
            r2 = reads[1]
        }

        tuple(meta, r1, r2)
    }


    //
    // MODULE: Assemble contigs with MEGAHIT
    // TODO: Check params
    MEGAHIT(ch_reads_host_depleted)

    ch_contigs = MEGAHIT.out.contigs
                            .multiMap { meta, contigs ->
                                meta: meta
                                contigs: contigs
                            }

    // ch_meta_long = ch_contigs.meta.collect()
    ch_contigs_long = ch_contigs.contigs.collect()

    //
    // MODULE: Concatenate contigs using a local module
    CONCATENATE_CONTIG_FILES(ch_contigs_long)



    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'viral_metagenomics_lite_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'viral_metagenomics_lite'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
