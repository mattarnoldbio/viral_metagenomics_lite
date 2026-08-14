#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mattarnoldbio/viral_metagenomics_lite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mattarnoldbio/viral_metagenomics_lite
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VIRAL_METAGENOMICS_LITE  } from './workflows/viral_metagenomics_lite'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_viral_metagenomics_lite_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_viral_metagenomics_lite_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_viral_metagenomics_lite_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    VIRAL_METAGENOMICS_LITE (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
        multiqc_report                 = VIRAL_METAGENOMICS_LITE.out.multiqc_report

        trimmed_and_deduplicated_reads = VIRAL_METAGENOMICS_LITE.out.trimmed_and_deduplicated_reads
        host_depleted_reads            = VIRAL_METAGENOMICS_LITE.out.host_depleted_reads

        separate_contigs               = VIRAL_METAGENOMICS_LITE.out.separate_contigs
        combined_contigs               = VIRAL_METAGENOMICS_LITE.out.combined_contigs
        all_virus_contigs              = VIRAL_METAGENOMICS_LITE.out.all_virus_contigs

        diamond_raw                    = VIRAL_METAGENOMICS_LITE.out.diamond_raw
        diamond_separate               = VIRAL_METAGENOMICS_LITE.out.diamond_separate
        diamond_combined               = VIRAL_METAGENOMICS_LITE.out.diamond_combined
        blastn_combined                = VIRAL_METAGENOMICS_LITE.out.blastn_combined
        blastn_separate_virus_hits     = VIRAL_METAGENOMICS_LITE.out.blastn_separate_virus_hits
        blastn_separate_non_virus_hits = VIRAL_METAGENOMICS_LITE.out.blastn_separate_non_virus_hits
        blastn_combined_virus_hits     = VIRAL_METAGENOMICS_LITE.out.blastn_combined_virus_hits
        blastn_combined_non_virus_hits = VIRAL_METAGENOMICS_LITE.out.blastn_combined_non_virus_hits

        final_filtered_virus_hits     = VIRAL_METAGENOMICS_LITE.out.final_filtered_virus_hits
        final_filtered_non_virus_hits = VIRAL_METAGENOMICS_LITE.out.final_filtered_non_virus_hits
    
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.multiqc_report
    )


    publish:
        multiqc_report = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.multiqc_report

        trimmed_and_deduplicated_reads = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.trimmed_and_deduplicated_reads
        host_depleted_reads            = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.host_depleted_reads

        separate_contigs  = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.separate_contigs
        combined_contigs  = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.combined_contigs
        all_virus_contigs = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.all_virus_contigs

        diamond_raw                     = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.diamond_raw
        diamond_separate                = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.diamond_separate
        diamond_combined                = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.diamond_combined
        blastn_combined                 = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.blastn_combined
        blastn_separate_virus_hits      = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.blastn_separate_virus_hits
        blastn_separate_non_virus_hits  = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.blastn_separate_non_virus_hits
        blastn_combined_virus_hits      = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.blastn_combined_virus_hits
        blastn_combined_non_virus_hits  = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.blastn_combined_non_virus_hits

        final_filtered_virus_hits     = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.final_filtered_virus_hits
        final_filtered_non_virus_hits = MATTARNOLDBIO_VIRAL_METAGENOMICS_LITE.out.final_filtered_non_virus_hits
}

output {
    multiqc_report {}

    trimmed_and_deduplicated_reads { path 'preprocessed_reads/trimmed_and_deduplicated_reads' }
    host_depleted_reads            { path 'preprocessed_reads/host_depleted_reads' }

    separate_contigs  { path 'contigs/separate_contigs' }
    combined_contigs  { path 'contigs/combined_contigs' }
    all_virus_contigs { path 'contigs/all_virus_contigs' }

    diamond_raw                     { path 'db_search_results' }
    diamond_separate                { path 'db_search_results' }
    diamond_combined                { path 'db_search_results' }
    blastn_combined                 { path 'db_search_results' }
    blastn_separate_virus_hits      { path 'db_search_results' }
    blastn_separate_non_virus_hits  { path 'db_search_results' }
    blastn_combined_virus_hits      { path 'db_search_results' }
    blastn_combined_non_virus_hits  { path 'db_search_results' }

    final_filtered_virus_hits     { path 'final_results' }
    final_filtered_non_virus_hits { path 'final_results' }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
