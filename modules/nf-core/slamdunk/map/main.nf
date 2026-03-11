process SLAMDUNK_MAP {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/slamdunk:0.4.3--py_0'
        : 'biocontainers/slamdunk:0.4.3--py_0'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${input}" == "${prefix}.bam") {
        error("Input and output names of bam files are the same, set prefix in module configuration to disambiguate!")
    }
    """
    slamdunk \\
        map \\
        -r ${fasta} \\
        -t ${task.cpus} \\
        -o outputs \\
        ${args} \\
        ${input}

    mv outputs/*.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slamdunk: \$(echo \$(slamdunk --version) | sed 's/^slamdunk //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${input}" == "${prefix}.bam") {
        error("Input and output names of bam files are the same, set prefix in module configuration to disambiguate!")
    }
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slamdunk: \$(echo \$(slamdunk --version) | sed 's/^slamdunk //')
    END_VERSIONS
    """
}
