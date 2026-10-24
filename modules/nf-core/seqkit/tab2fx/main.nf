// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SEQKIT_TAB2FX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/297e985dd395d87f9089b314dbe685ddce59fa0770c4d1f5e1009bc7b6ec7a8d/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835'}"

    input:
    tuple val(meta), path(text)
    val out_ext

    output:
    tuple val(meta), path("*.f*"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "fa.zst"
    """
    seqkit \\
        tab2fx \\
        ${args} \\
        --threads ${task.cpus} \\
        ${text} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "fa.zst"
    """
    echo ${args}

    touch ${prefix}.${suffix}
    """
}
