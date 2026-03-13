process CNVPYTOR_PARTITION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bbb6343edff4191cb1f445b2aac028d1f805ed5a7d50799513c82531bcfdede5/data'
        : 'community.wave.seqera.io/library/cnvpytor_make:a8fdcebe82041114'}"

    input:
    tuple val(meta), path(pytor)
    val bin_sizes

    output:
    tuple val(meta), path("${prefix}.pytor"), emit: pytor
    // cnvpytor version is hardcoded due to this error when calling cnvpytor --version
    // > cnvpytor --version
    // 2026-03-13 16:14:08,088 - cnvpytor - ERROR - Some reference genome resource files are missing.
    // Run 'cnvpytor -download' as same user who has installed cnvpytor.
    tuple val("${task.process}"), val('cnvpytor'), val('1.3.2'), emit: versions_cnvpytor, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bins = bin_sizes ?: '1000'
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(pwd)/.mplconfig

    cnvpytor \\
        -root ${pytor} \\
        -partition ${bins} \\
        ${args}

    cp ${pytor} ${prefix}.pytor
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor
    """
}
