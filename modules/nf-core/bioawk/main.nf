process BIOAWK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6' :
        'quay.io/biocontainers/bioawk:1.0--h5bf99c6_6' }"

    input:
    tuple val(meta), path(input)
    path(program_file)
    val(disable_redirect_output)
    val output_file_extension

    output:
    tuple val(meta), path("*.${output_file_extension}"), emit: output
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.
    tuple val("${task.process}"), val('bioawk'), val("1.0"), emit: versions_bioawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    program         = program_file ? "-f ${program_file}" : "${args2}"
    def prefix      = task.ext.prefix ?: "${meta.id}"
    output_cmd      = output_file_extension.endsWith("gz") ? "| gzip > ${prefix}.${output_file_extension}" : "> ${prefix}.${output_file_extension}"
    output          = disable_redirect_output ? "" : output_cmd

    if ("${input}" == "${prefix}.${output_file_extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate."

    """
    bioawk \\
        ${args} \\
        ${program} \\
        ${input} \\
        ${output}
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def create_cmd = output_file_extension.endsWith("gz") ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${output_file_extension}
    """
}
