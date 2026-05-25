process HAMRONIZATION_AMRFINDERPLUS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(report)
    val format
    val software_version
    val reference_db_version

    output:
    tuple val(meta), path("*.json"), optional: true, emit: json
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val("${task.process}"), val('hamronization'), eval("hamronize --version 2>&1 | sed 's/hamronize //'"), topic: versions, emit: versions_hamronization

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hamronize \\
        amrfinderplus \\
        ${report} \\
        ${args} \\
        --format ${format} \\
        --analysis_software_version ${software_version} \\
        --reference_database_version ${reference_db_version} \\
        --input_file_name ${prefix} \\
        > ${prefix}.${format}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${format}
    """
}
