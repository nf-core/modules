// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ICOUNTMINI_PEAKS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e0ecd2ae293566c4a35a7e022f42220c41449d02a40c690f3132955088ce2ce0/data':
        'community.wave.seqera.io/library/icount-mini_gzip:00407b12cb167b11' }"

    input:
    tuple val(meta), path(bed), path(sigxls)

    output:
    tuple val(meta), path("*.peaks.bed.gz"), emit: peaks
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini peaks \\
        ${bed} \\
        ${sigxls} \\
        ${prefix}.peaks.bed.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.peaks.bed.gz
    """
}
