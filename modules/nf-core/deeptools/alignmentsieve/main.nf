process DEEPTOOLS_ALIGNMENTSIEVE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:41defd13a6f2ce014549fcc05d0b051f655777f9-0':
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:41defd13a6f2ce014549fcc05d0b051f655777f9-0' }"

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*_as.bam") , emit: bam
    path  "versions.yml"              , emit: versions
    path  "*_log.txt"                 , emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    alignmentSieve \\
        $args \\
        -b $input \\
        -o ${prefix}_as.bam \\
        --filterMetrics ${prefix}_log.txt \\
        --numberOfProcessors $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(alignmentSieve --version | sed -e "s/alignmentSieve //g")
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_as.bam
    touch ${prefix}_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(alignmentSieve --version | sed -e "s/alignmentSieve //g")
    END_VERSIONS
    """
}
