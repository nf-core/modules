process QSV_CAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/864fcf1e6a250b588eaabd4df01334e44bc528fd29049c4fb24f7c56c8c33695/data'
        : 'community.wave.seqera.io/library/qsv:5.1.0--9a6a0c23d3b279b5'}"

    input:
    tuple val(meta), path(csv, name: 'inputs/in*/*')
    val mode
    val out_format
    val skip_input_format_check

    output:
    tuple val(meta), path("${prefix}.${out_format}"), emit: csv
    tuple val("${task.process}"), val('qsv'), eval("qsv --version | cut -d' ' -f2 | cut -d'-' -f1"), topic: versions, emit: versions_qsv

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def skip_format_check_cmd = skip_input_format_check ? 'export QSV_SKIP_FORMAT_CHECK=1' : ''
    """
    ${skip_format_check_cmd}

    qsv \\
        cat \\
        ${mode} \\
        ${args} \\
        -o ${prefix}.${out_format} \\
        ${csv}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.${out_format}
    """
}
