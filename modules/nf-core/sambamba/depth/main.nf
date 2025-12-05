process SAMBAMBA_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--he614052_3':
        'biocontainers/sambamba:1.0.1--he614052_3' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(bed)
    val(mode)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!['region','window','base'].contains(mode)) {
        error "Mode needs to be one of: region, window, base"
    }
    def bed_arg = bed ? "--regions ${bed}" : ''

    """
    sambamba \\
        depth \\
        $mode \\
        $bed_arg \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}.bed \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}
