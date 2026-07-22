process SPLITBLASTRESULT {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    tuple val(meta), val(blast_results)

    output:
    tuple val(meta), path("*.txt")

    script:
    def sample = meta.id
    def db = meta.db
    def file_contents = blast_results.collect{ row -> row.join("\t")}.join("\n")
    """
    touch ${sample}.${db}.txt
    echo -e "${file_contents}" >> ${sample}.${db}.txt

    """
}