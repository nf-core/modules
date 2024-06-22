process MINDAGAP_DUPLICATEFINDER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mindagap=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mindagap:0.0.2--pyhdfd78af_1':
        'biocontainers/mindagap:0.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(spot_table)

    output:
    tuple val(meta), path("*markedDups.txt"), emit: marked_dups_spots
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    duplicate_finder.py \\
        $spot_table \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mindagap: \$(mindagap.py test -v)
    END_VERSIONS
    """
}
