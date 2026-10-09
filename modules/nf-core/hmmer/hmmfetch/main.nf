// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HMMER_HMMFETCH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4'
        : 'quay.io/biocontainers/hmmer:3.4--hb6cb901_4'}"

    input:
    tuple val(meta), path(hmm)
    val key
    path keyfile
    path index

    output:
    tuple val(meta), path("*.hmm"), emit: hmm, optional: true
    tuple val(meta), path("*.ssi"), emit: index, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keyarg = key ?: ''
    def kfopt = keyfile ? '-f' : ''
    def index_opt = !key && !keyfile ? '--index' : ''
    def outfile = !key && !keyfile ? '' : "> ${prefix}.hmm"

    """
    hmmfetch \\
        ${kfopt} \\
        ${index_opt} \\
        ${args} \\
        ${hmm} \\
        ${keyarg} \\
        ${keyfile} \\
        ${outfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.hmm
    """
}
