process ASSIGN_TAXONOMY {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    tuple val(meta), val(blast_results)
    val which_db
    tuple val(rankedlineage), val(fullnamelineage)

    output:
    tuple val(meta), path("*.tsv")

    script:
    """
    AssignTaxonomy.sh $blast_results $which_db $rankedlineage $fullnamelineage 

    """
}