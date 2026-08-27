process ASSIGN_TAXONOMY {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.4.1--8dd22ef9018762ae' :
        'community.wave.seqera.io/library/gawk:5.4.1--8dd22ef9018762ae' }"

    input:
    tuple val(meta), path(blast_results)
    val which_db
    tuple path(rankedlineage), path(fullnamelineage)

    output:
    tuple val(meta), path("*.tsv")

    script:
    """
    AssignTaxonomy.sh $blast_results $which_db $rankedlineage $fullnamelineage 

    """
}