/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                   } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE                                                } from '../modules/nf-core/trimgalore/main'
include { PRINSEQPLUSPLUS                                           } from '../modules/nf-core/prinseqplusplus/main'
include { BOWTIE2_BUILD                                             } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                                             } from '../modules/nf-core/bowtie2/align/main'
include { MEGAHIT                                                   } from '../modules/nf-core/megahit/main'
include { CONCATENATECONTIGFILES                                    } from '../modules/local/filewrangling/concatenatecontigfiles/main' 
include { SPLITBLASTRESULT as SPLITBLASTRESULT_BLASTX               } from '../modules/local/filewrangling/splitblastresult/main' 
include { SPLITBLASTRESULT as SPLITBLASTRESULT_BLASTN               } from '../modules/local/filewrangling/splitblastresult/main' 
include { DIAMOND_BLASTX                                            } from '../modules/nf-core/diamond/blastx/main'
include { KRONA_KTIMPORTBLAST as KRONA_KTIMPORTBLAST_BLASTX         } from '../modules/local/kronatools/importblast/main'
include { KRONA_KTIMPORTBLAST as KRONA_KTIMPORTBLAST_BLASTN         } from '../modules/local/kronatools/importblast/main'
include { EXTRACT_KRONA_TAXONOMY as EXTRACT_KRONA_TAXONOMY_BLASTX   } from '../modules/local/extractkronataxonomy/main'
include { EXTRACT_KRONA_TAXONOMY as EXTRACT_KRONA_TAXONOMY_BLASTN   } from '../modules/local/extractkronataxonomy/main'
include { CONCATENATECSVS as CONCATENATECSVS_BLASTX                 } from '../modules/local/filewrangling/concatenatecsvs/main'
include { CONCATENATECSVS as CONCATENATECSVS_BLASTN                 } from '../modules/local/filewrangling/concatenatecsvs/main'
include { EXTRACTVIRUSCONTIGS as EXTRACTVIRUSCONTIGS_BLASTX         } from '../modules/local/filewrangling/extractviruscontigs/main'
include { EXTRACTVIRUSCONTIGS as EXTRACTVIRUSCONTIGS_BLASTN         } from '../modules/local/filewrangling/extractviruscontigs/main'
include { CONCATENATECONTIGFILES as  CONCATENATECONTIGFILES_BLASTX  } from '../modules/local/filewrangling/concatenatecontigfiles/main' 
include { BLAST_BLASTN                                              } from '../modules/nf-core/blast/blastn/main'
include { paramsSummaryMap                                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                    } from '../subworkflows/local/utils_nfcore_viral_metagenomics_lite_pipeline'

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
    ch_trimmed = TRIMGALORE.out.reads

    //
    // MODULE: Remove PCR duplicates and low quality reads with PRINSEQ++
    // TODO: Check params
    PRINSEQPLUSPLUS(ch_trimmed)
    ch_dedup = PRINSEQPLUSPLUS.out.good_reads
    
    ch_dedup
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

    MEGAHIT.out.contigs
               .multiMap { meta, contigs ->
                    meta: meta
                    contigs: [ [id: meta.id], contigs ]
                    contigs_: contigs
                }
                .set { ch_contigs }


    // ch_meta_long = ch_contigs.meta.collect()
    ch_contigs_long = ch_contigs.contigs_.collect()

    //
    // MODULE: Concatenate contigs using a local module
    CONCATENATECONTIGFILES(ch_contigs_long, "")

    ch_all_contigs = CONCATENATECONTIGFILES.out.map{ all_contigs -> 
        def meta = [id: "all_contigs"]
        tuple (meta, all_contigs)
    }
    ch_diamond_db = channel.fromPath(params.diamond_db, checkIfExists: true)
                            .map{
                                db_path -> 
                                def meta = [id: "blast_db"]

                                tuple(meta, db_path)
                            }


    //
    // MODULE: Search contigs against reference DB with Diamond BLASTX
    DIAMOND_BLASTX(ch_all_contigs, ch_diamond_db ,'txt', '')


    // ch_krona_db = channel.fromPath(params.krona_db, checkIfExists: true)
    //                     .map{
    //                         db_path -> 
    //                         def meta = [id: "blast_db"]

    //                         tuple(meta, db_path)
    //                     }

    ch_blast_results = DIAMOND_BLASTX.out.txt
                        .flatMap{
                            _meta, results ->
                                results.splitCsv(sep: "\t")
                                       .collect { row ->
                                             def sample_id = row[0].split('\\|')[0]   // parse "sample" out of "sample|contig"
                                             tuple(sample_id, row)
                                            }
                        }
                        .groupTuple()
                        .map { sample_id, blast_results ->
                                  def meta = [id: sample_id, db: "diamond"]
                               tuple(meta, blast_results)
                                }
    
    // MODULE: Split BLAST results into individual files for each sample
    SPLITBLASTRESULT_BLASTX(ch_blast_results)
    ch_split_blastx_results = SPLITBLASTRESULT_BLASTX.out
    // ch_split_blastx_results.view()

    //
    // MODULE: Make Krona plots to show taxonomic breakdown of each sample's BLAST results
    KRONA_KTIMPORTBLAST_BLASTX(ch_split_blastx_results, params.krona_db)

    ch_contigs_taxonomy = ch_contigs.contigs
        .map { meta, contigs -> tuple(meta.id, contigs) }
        .join(
            KRONA_KTIMPORTBLAST_BLASTX.out.html.map { meta, html -> tuple(meta.id, meta, html) }
        )
        .map { _id, contigs, meta, html -> tuple(meta, contigs, html) }
    //if legacy :
    //TODO: Add a parameter to choose between legacy and new taxonomy extraction
    //TODO: Add a parameter to choose between database and score filter for Krona taxonomy extraction

    // MODULE: Extract virus hits from Krona taxonomy results
    EXTRACT_KRONA_TAXONOMY_BLASTX(ch_contigs_taxonomy, "diamond", 10)

    EXTRACT_KRONA_TAXONOMY_BLASTX.out.virus_hits
        .set{ ch_virus_hits }

    EXTRACT_KRONA_TAXONOMY_BLASTX.out.virus_hits
                        .map { _meta, _contigs, virus_hits ->
                            virus_hits
                        }
                        .collect()
                        .set { ch_all_virus_hits }
    
    // MODULE: Concatenate .csv files from EXTRACT_KRONA_TAXONOMY_BLASTX to create a single .csv file for virus hits
    CONCATENATECSVS_BLASTX(ch_all_virus_hits, "diamond", "all")

    // MODULE: Extract virus contigs from the original contigs file based on the virus hits .csv file
    EXTRACTVIRUSCONTIGS_BLASTX(ch_virus_hits)
    
    EXTRACTVIRUSCONTIGS_BLASTX.out.virus_contigs
        .set{ ch_all_virus_contigs }

    ch_all_virus_contigs
        .map({ _meta, virus_contigs ->
            virus_contigs
        })
        .collect()
        .set { ch_all_virus_contigs_flat }

    // MODULE: Concatenate virus contigs from EXTRACTVIRUSCONTIGS_BLASTX to create a single .fa.gz file for virus contigs
    CONCATENATECONTIGFILES_BLASTX(ch_all_virus_contigs_flat, "_diamond_virus_hits")
        .map{ all_virus_contigs ->
            def meta = [id: "all_virus_contigs"]
            tuple(meta, all_virus_contigs)

        }        
        .set { ch_all_virus_contigs_fasta }


    ch_blast_db = channel.fromPath(params.blast_db, checkIfExists: true)
                        .map{
                            db_path -> 
                            def meta = [id: "blast_db"]

                            tuple(meta, db_path)
                        }

    // MODULE: Run BLASTN search on the virus hits against nt    
    BLAST_BLASTN(ch_all_virus_contigs_fasta, ch_blast_db, [], [], false)
        
    ch_blastn_results = BLAST_BLASTN.out.txt
                        .flatMap{
                            _meta, results ->
                                results.splitCsv(sep: "\t")
                                       .collect { row ->
                                             def sample_id = row[0].split('\\|')[0]   // parse "sample" out of "sample|contig"
                                             tuple(sample_id, row)
                                            }
                        }
                        .groupTuple()
                        .map { sample_id, blast_results ->
                                  def meta = [id: sample_id, db: "blastn"]
                               tuple(meta, blast_results)
                                }

    SPLITBLASTRESULT_BLASTN(ch_blastn_results)

    SPLITBLASTRESULT_BLASTN.out
        .map{ _meta, blast_results ->
            def meta = [id: _meta.id, db: "blastn"]
            tuple(meta, blast_results)
        }
        .set { ch_split_blastn_results }


    KRONA_KTIMPORTBLAST_BLASTN(ch_split_blastn_results, params.krona_db)


    ch_blast_contigs_taxonomy = ch_all_virus_contigs
        .map { meta, contigs -> tuple(meta.id, contigs) }
        .join(
            KRONA_KTIMPORTBLAST_BLASTN.out.html.map { meta, html -> tuple(meta.id, meta, html) }
        )
        .map { _id, contigs, meta, html -> tuple(meta, contigs, html) }


    //if legacy :
    //TODO: Add a parameter to choose between legacy and new taxonomy extraction
    //TODO: Add a parameter to choose between database and score filter for Krona taxonomy extraction

    // MODULE: Extract virus hits from Krona taxonomy results
    EXTRACT_KRONA_TAXONOMY_BLASTN(ch_blast_contigs_taxonomy, "blastn", 10)

    EXTRACT_KRONA_TAXONOMY_BLASTN.out.non_virus_hits
        .map{ _meta, _contigs, non_virus_hits ->
            non_virus_hits
        }
        .collect()
        .set{ ch_non_virus_hits }

    ch_non_virus_hits.view()    

    // MODULE: Concatenate .csv files from EXTRACT_KRONA_TAXONOMY_BLASTN to create a single .csv file for non-virus hits
    CONCATENATECSVS_BLASTN(ch_non_virus_hits, "blastn", "non")

    

    //TODO: classify the BLASTN results

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
