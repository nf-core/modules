process MOTUS_PREPLONG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.gz"), emit: out
    // WARN: Version information not provided by tool on CLI.  Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('motus'), val("3.1.0"), topic: versions, emit: versions_motus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def refdb = db ? "-db ${db}" : ""

    """
    motus \\
        prep_long \\
        ${args} \\
        -i ${reads} \\
        ${refdb} \\
        -t ${task.cpus} \\
        -o ${prefix}.gz \\
        2>| >(tee ${prefix}.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.gz
    """
}
