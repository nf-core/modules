process SKA_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska:1.0--h077b44d_6':
        'biocontainers/ska:1.0--h077b44d_6' }"

    input:
    tuple val(meta), path(sketch_files), path(sketch_list)

    output:
    tuple val(meta), path("*distances.tsv"), emit: distances    , optional: true
    tuple val(meta), path("*clusters.tsv") , emit: cluster_list , optional: true
    tuple val(meta), path("*cluster*.txt") , emit: cluster_files, optional: true
    tuple val(meta), path("*.dot")         , emit: dot          , optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_list = sketch_list ? "-f ${sketch_list}" : ''
    """
    ska \\
        distance \\
        $args \\
        $arg_list \\
        -o ${prefix} \\
        $sketch_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version | grep Version |& sed '1!d ; s/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.distances.tsv
    touch ${prefix}.clusters.tsv
    touch ${prefix}.cluster.1.txt
    touch ${prefix}.dot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version | grep Version |& sed '1!d ; s/Version: //')
    END_VERSIONS
    """
}
