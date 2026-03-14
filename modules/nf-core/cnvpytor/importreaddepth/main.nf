process CNVPYTOR_IMPORTREADDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bb/bbb6343edff4191cb1f445b2aac028d1f805ed5a7d50799513c82531bcfdede5/data':
        'community.wave.seqera.io/library/cnvpytor_make:a8fdcebe82041114' }"

    input:
    tuple val(meta), path(input_file), path(index)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.pytor")	, emit: pytor
    tuple val("${task.process}"), val('cnvpytor'), val('1.3.2'), emit: versions_cnvpytor, topic: versions
    // cnvpytor version is hardcoded due to this error when calling cnvpytor --version
    // > cnvpytor --version
    // 2026-03-13 16:14:08,088 - cnvpytor - ERROR - Some reference genome resource files are missing.
    // Run 'cnvpytor -download' as same user who has installed cnvpytor.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "-T ${fasta}" : ''
    """
    cnvpytor \\
        -root ${prefix}.pytor \\
        -rd $input_file \\
        $args \\
        $reference
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor
    """
}
