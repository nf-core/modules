// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AGAT_SPKEEPLONGESTISOFORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdf86a4db56d388b2cad1599b640e1e33ab964dcdb74a8147981d48a429bc7d5/data'
:         'community.wave.seqera.io/library/agat:1.6.1--1f8d94d325805bb4' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ""
    output           = "${prefix}.longest.gff"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff ${gxf} \\
        --cpu ${task.cpus} \\
        ${config_param} \\
        --out ${output} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output     = "${prefix}.longest.gff"
    """
    touch ${output}
    """
}
