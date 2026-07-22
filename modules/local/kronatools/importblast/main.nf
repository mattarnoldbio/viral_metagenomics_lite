process KRONA_KTIMPORTBLAST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1' :
        'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(blast_output)
    path db

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.db}"
    """
    ktImportBLAST \\
        $args \\
        -o ${prefix}.html \\
        $blast_output \\
        -tax $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(echo \$(ktImportBLAST 2>&1) | grep -Po '(?<=KronaTools )[0-9.]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.db}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(echo \$(ktImportBLAST 2>&1) | grep -Po '(?<=KronaTools )[0-9.]+')
    END_VERSIONS
    """
}