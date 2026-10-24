// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process WGSIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e41ff04a2e8632e62ecbce8ac2e3e06f6738c815b1a4f38acd937b133c587b92/data'
:         'community.wave.seqera.io/library/wgsim:1.0--2fb08e23ce8d62af' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    tuple val("${task.process}"), val('wgsim'), eval("wgsim 2>&1 | sed -n 's/Version: //p'"), topic: versions, emit: versions_wgsim

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    wgsim \\
        ${args} \\
        ${fasta} \\
        ${prefix}_R1.fastq \\
        ${prefix}_R2.fastq
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq
    touch ${prefix}_R2.fastq
    """
}
