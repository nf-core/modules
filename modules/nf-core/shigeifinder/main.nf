process SHIGEIFINDER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigeifinder:1.3.5--pyhdfd78af_0':
        'quay.io/biocontainers/shigeifinder:1.3.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('shigapass'), val("1.3.5"), topic: versions, emit: versions_shigapass

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shigeifinder \\
        ${args} \\
        --output ${prefix}.tsv \\
        -t ${task.cpus} \\
        -i ${seqs}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
