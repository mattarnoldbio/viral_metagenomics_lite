process FILTERHITDATA {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/beautifulsoup4_pandas_pigz:c51f352c49b97995' :
        'community.wave.seqera.io/library/beautifulsoup4_pandas_pigz:c51f352c49b97995' }"

    input:
    path(diamond_virus_hits)
    path(blastn_virus_hits)
    path(non_virus_hits)
    tuple val(meta), path(all_virus_contigs_fasta)
    val(min_contig_length)

    output:
    path("finalvertebrate_viral_hits.csv"), optional: true,  emit: vertebrate_virus_hits
    path("finalnon_vertebrate_viral_hits.csv"), optional: true,  emit: non_vertebrate_virus_hits

    script:
    """
    contigs_file=${all_virus_contigs_fasta}
    uncompressed_contigs=\${contigs_file%".gz"}
    pigz -dc \$contigs_file > \$uncompressed_contigs

    FilterDiamondHits.py -d ./ --diamond_cutoff 0 --blast_cutoff 0 --min_contig_length $min_contig_length --filter_vertebrate

    touch finalvertebrate_viral_hits.csv
    touch finalnon_vertebrate_viral_hits.csv

    """
}