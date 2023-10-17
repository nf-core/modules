process WINDOWMASKER_MKCOUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::blast=2.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.0--h7d5a4b4_1':
        'biocontainers/blast:2.14.0--h7d5a4b4_1' }"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.txt")  , emit: counts
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    def memory = 3072
    if (!task.memory) {
        log.info '[WINDOWMASKER: MK_COUNTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        memory = (task.memory.toMega()).intValue()
    }

    """
    windowmasker -mk_counts \\
        $args \\
        -mem ${memory} \\
        -in ${ref} \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """
}
