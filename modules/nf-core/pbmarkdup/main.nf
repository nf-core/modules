process PBMARKDUP {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmarkdup:1.2.0--h9ee0642_0' :
        'biocontainers/pbmarkdup:1.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.{bam,f*a,/.*f.*\\.gz/}") , emit: markduped
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = input.getExtension()

    """
    pbmarkdup \\
        -j ${task.cpus} \\
        $input \\
        ${prefix}.${suffix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmarkdup: \$(echo \$(pbmarkdup --version 2>&1) | awk 'BEFORE{FS=" "}{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = input.getExtension()
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmarkdup: \$(echo \$(pbmarkdup --version 2>&1) | awk 'BEFORE{FS=" "}{print \$2}')
    END_VERSIONS
    """
}
