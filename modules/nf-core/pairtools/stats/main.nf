process PAIRTOOLS_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pairtools=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0' :
        'quay.io/biocontainers/pairtools:1.0.2--py39h2a9f597_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.stats.tsv")   , emit: stat
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
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
