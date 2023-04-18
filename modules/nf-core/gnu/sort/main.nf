process GNU_SORT {
    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), file( "*.sorted" )   , emit: sorted
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def output      = task.ext.prefix   ?: "${meta.id}.txt"
    output_file     = "${output}.sorted"
    def VERSION     = "9.1"             // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    sort ${args} ${file} >  ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args     ?: ''
    def output      = task.ext.prefix   ?: "${meta.id}"
    output_file     = "${output}.sorted"
    def VERSION     = "9.1"
    """
    cat ${file} | sort ${args} > ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
