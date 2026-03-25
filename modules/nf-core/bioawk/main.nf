process BIOAWK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6' :
        'biocontainers/bioawk:1.0--h5bf99c6_6' }"

    input:
    tuple val(meta), path(input)
    val suffix
    val zip_bool

    output:
    tuple val(meta), path("${file_output}"),    optional: true, emit: output
    tuple val(meta), path("*.gz"),              optional: true, emit: gz_output
    tuple val("${task.process}"), val('bioawk'), val("1.0"), emit: versions_bioawk, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update version string above when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    file_output     = "${prefix}.${suffix}"
    if ("${input}" == "${file_output}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate."

    """
    bioawk \
        $args \
        $input \
        > ${file_output}

    if [ "${zip_bool}" = "true" ]; then
        gzip "${file_output}"
    fi

    """

    stub:
    """
    touch ${file_output}
    echo "" > ${file_output}.gz
    """
}
