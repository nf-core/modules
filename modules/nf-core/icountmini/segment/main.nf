// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ICOUNTMINI_SEGMENT {
    tag "$gtf"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/0905aab0a42e3ebef7c4300c7d8fdf054c63e805cf2a1cae45027658a5c4a88a/data':
        'community.wave.seqera.io/library/icount-mini:4.0.2--4fff630fa45c91ab' }"

    input:
    tuple val(meta), path(gtf)
    path fai

    output:
    tuple val(meta), path("*_seg.gtf")       ,  emit: gtf
    tuple val(meta), path("*_regions.gtf.gz"),  emit: regions
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"
    """
    iCount-Mini segment \\
        $args \\
        $gtf \\
        ${prefix}_seg.gtf \\
        $fai

    mv regions.gtf.gz ${prefix}_regions.gtf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.simpleName}"
    """
    touch ${prefix}_seg.gtf
    echo "" | gzip > ${prefix}_regions.gtf.gz
    """
}
