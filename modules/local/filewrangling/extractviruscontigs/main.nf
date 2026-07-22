process EXTRACTVIRUSCONTIGS {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    tuple val(meta), path(contigs), path(virus_hits_csv)

    output:
    tuple val(meta), path("*_diamond_virus_hits.fa.gz") ,emit: virus_contigs

    script:
    """

    touch ${meta.id}_diamond_virus_hits.fa

    virus_hits=\$(cat $virus_hits_csv | cut -f 2 -d "," | tail -n+2 ) 
    
    contigs_file=$contigs
    uncompressed_contigs="\${contigs_file%.gz}"
    pigz -dc \$contigs_file > \$uncompressed_contigs

    for hit in \$virus_hits; do # For each hit
        contig=\$(grep "\${hit} " \$uncompressed_contigs -A 1) 
        echo "\$contig" >> ${meta.id}_diamond_virus_hits.fa 
    done

    pigz ${meta.id}_diamond_virus_hits.fa

    """
}