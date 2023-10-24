process GOLEFT_INDEXSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::goleft=0.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goleft:0.2.4--h9ee0642_1':
        'biocontainers/goleft:0.2.4--h9ee0642_1' }"

    input:
    tuple val(meta) , path(bai)
    tuple val(meta2), path(fai)
    val split

    output:
    tuple val(meta), path("*.bed") , emit: bed
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    goleft \\
        indexsplit \\
        $args \\
        -n ${split} \\
        --fai ${fai} \\
        ${bai} \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: \$(goleft --version 2>&1 | head -n 1 | sed 's/^.*goleft Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: \$(goleft --version 2>&1 | head -n 1 | sed 's/^.*goleft Version: //')
    END_VERSIONS
    """
}
