process PAIRTOOLS_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pairtools=0.3.0 numpy=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5' :
        'quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5' }"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path("*.stats.tsv")   , emit: stats
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools \\
        stats \\
        $args \\
        --nproc-in $task.cpus \\
        --nproc-out $task.cpus \\
        -o ${prefix}.stats.tsv \\
        $pairs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
