
process PNEUMOCAT {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pneumocat:1.2.1--0':
        'quay.io/biocontainers/pneumocat:1.2.1--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.xml"), emit: xml
    tuple val(meta), path("*.txt"), emit: txt
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('pneumocat'), val("1.2.1"), topic: versions, emit: versions_pneumocat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        ${args} \\
        --threads ${task.cpus} \\
        --output_dir .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.results.xml
    touch ${prefix}.txt
    """
}
