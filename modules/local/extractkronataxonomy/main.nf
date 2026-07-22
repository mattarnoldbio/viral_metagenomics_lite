process EXTRACT_KRONA_TAXONOMY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beautifulsoup4_pandas_pigz:c51f352c49b97995' :
        'community.wave.seqera.io/library/beautifulsoup4_pandas_pigz:c51f352c49b97995' }"

    input:
    tuple val(meta), path(contigs), path(taxonomy_file)
    val which_db 
    val score_filter

    output:
    tuple val(meta), path(contigs) , path("*all_virus_hits.csv"), emit: virus_hits
    tuple val(meta), path(contigs) , path("*_non_virus_hits.csv"), optional: true, emit: non_virus_hits
    tuple val(meta), path(contigs) , path("*_no_virus_hits.csv"), optional: true, emit: no_virus_hits

    script:
    """
    contigs_file=$contigs
    uncompressed_contigs="\${contigs_file%.gz}"
    pigz -dc \$contigs_file > \$uncompressed_contigs
    ParseKrona.py -k $taxonomy_file -c \$uncompressed_contigs -o ./ -a $meta.id -w $which_db -s $score_filter --contig_mode
    """
}