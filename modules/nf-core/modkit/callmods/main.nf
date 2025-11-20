process MODKIT_CALLMODS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ont-modkit:0.5.0--hcdda2d0_2'
        : 'biocontainers/ont-modkit:0.5.0--hcdda2d0_2'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    modkit \\
        call-mods \\
        ${args} \\
        ${bam} \\
        ${prefix}.bam \\
        --threads ${task.cpus} \\
        --log-filepath ./${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/modkit //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.bam
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/modkit //' )
    END_VERSIONS
    """
}
