process FASTQE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastqe:0.5.2--pyhdfd78af_0'
        : 'biocontainers/fastqe:0.5.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastqe'), val("0.5.2"), emit: versions_fastqe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqe \\
        ${args} \\
        ${fastq} \\
        --output ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
