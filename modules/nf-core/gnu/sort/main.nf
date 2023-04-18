process GNU_SORT {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), file( "${output_file}" )   , emit: sorted
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: "${input.extension}"
    output_file     = "${prefix}.${suffix}"
    def VERSION     = "9.1"             // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    sort ${args} ${input} > ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: "${input.extension}"
    output_file     = "${prefix}.${suffix}"
    def VERSION     = "9.1"
    """
    sort ${args} ${input} > ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
