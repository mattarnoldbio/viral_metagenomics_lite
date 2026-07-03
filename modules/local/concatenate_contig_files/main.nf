process CONCATENATE_CONTIG_FILES {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    path(contigs)

    output:
    path("all_contigs.fa.gz"), emit: contigs

    script:
    """
    contigs=(${contigs})

    touch all_contigs.fa

    for contigs_file in "\${contigs[@]}" ; do
        uncompressed_contigs=\${contigs_file%".gz"}
        accession=\${contigs_file%".contigs.fa.gz"}
        pigz -dc \$contigs_file > \$uncompressed_contigs
        cat \$uncompressed_contigs | sed "s/^>/>\${accession}|/" >> all_contigs.fa
    done

    pigz all_contigs.fa

    """
}