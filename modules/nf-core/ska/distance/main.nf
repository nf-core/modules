process SKA_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska:1.0--h077b44d_6':
        'biocontainers/ska:1.0--h077b44d_6' }"

    input:
    tuple val(meta), path(sketch_files, arity: '0..*'), path(sketch_list)

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
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def output_dist  = task.ext.args =~ "-d" ? "" : "touch ${prefix}.distances.tsv"
    def output_clust = task.ext.args =~ "-c" ? "" : "touch ${prefix}.clusters.tsv"
    // this is not a complete criterion for this output but it is good enough
    def output_dot   = task.ext.args =~ "-S" || sketch_files.size > 1 ? "touch ${prefix}.dot" : ""
    """
    $output_dist
    $output_clust
    # this is not how this works but it's the best we can do without knowing the input content
    for i in {1..${sketch_files.size}}
    do
        touch ${prefix}.cluster\${i}.txt
    done
    $output_dot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version | grep Version |& sed '1!d ; s/Version: //')
    END_VERSIONS
    """
}
