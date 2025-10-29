process CMAPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmaple:1.1.0--h503566f_1':
        'biocontainers/cmaple:1.1.0--h503566f_1' }"

    input:
    tuple val(meta), path(aln), path(newick)

    output:
    tuple val(meta), path("*.treefile"), emit: treefile
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def tree_arg = newick ? "-t ${newick}" : ""
    """
    cmaple-aa \\
        $args \\
        -nt $task.cpus \\
        --prefix ${prefix} \\
        ${tree_arg} \\
        -aln $aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmaple: \$(cmaple --help | grep -m1 'CMAPLE version' | sed -E 's/.*version ([0-9.]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.treefile
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmaple: \$(cmaple --help | grep -m1 'CMAPLE version' | sed -E 's/.*version ([0-9.]+).*/\\1/')
    END_VERSIONS
    """
}
