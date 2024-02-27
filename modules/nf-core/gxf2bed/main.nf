process GXF2BED {
    tag '$gxf.baseName'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gxf2bed:0.2.0--h4ac6f70_0':
        'biocontainers/gxf2bed:0.2.0--h4ac6f70_0 ' }"

    input:
    path gxf

    output:
    path "*.bed"       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${gxf.baseName}"
    def args = task.ext.args ?: ''
    """
    gxf2bed \\
        $args \\
        --input $gxf \\
        --threads $task.cpus \\
        --output ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$(gxf2bed --version |& sed '1!d ; s/gxf2bed //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gxf.baseName}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$(gxf2bed --version |& sed '1!d ; s/gxf2bed //')
    END_VERSIONS
    """
}
