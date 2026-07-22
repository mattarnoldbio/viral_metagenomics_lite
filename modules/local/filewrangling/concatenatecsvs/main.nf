process CONCATENATECSVS {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    conda "bioconda::pigz=2.3.4"
    input:
    path(csv_files)

    output:
    path("all_virus_hits.csv"), emit: all_virus_hits_csv

    script:
    """
    csvs=(${csv_files})

    touch all_virus_hits.csv

    for csv in "\${csvs[@]}" ; do
        tail -n+2 \${csv} >> all_virus_hits.csv
    done

    """
}