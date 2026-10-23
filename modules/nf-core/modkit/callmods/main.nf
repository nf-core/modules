// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MODKIT_CALLMODS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65e515e445d31be38188af5fe636e71af18abc601afacf80e9021a4752da6332/data'
:         'community.wave.seqera.io/library/ont-modkit:0.6.4--0c623277f8775a30' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    modkit \\
        call-mods \\
        ${args} \\
        ${bam} \\
        ${prefix}.bam \\
        --threads ${task.cpus} \\
        --log-filepath ./${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.bam
    touch ${prefix}.log
    """
}
