process CONCATENATECONTIGFILES {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    path(contigs)
    val(filename)

    output:
    path("all_contigs*.fa.gz"), emit: all_contigs

    script:
    """
    contigs=(${contigs})

    touch all_contigs${filename}.fa

    for contigs_file in "\${contigs[@]}" ; do
        uncompressed_contigs=\${contigs_file%".gz"}
        accession="\${contigs_file%%[._]*}"
        pigz -dc \$contigs_file > \$uncompressed_contigs
        cat \$uncompressed_contigs | sed "s/^>/>\${accession}|/" >> all_contigs${filename}.fa
    done

    pigz all_contigs${filename}.fa

    """
}