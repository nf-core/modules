// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CDHIT_CDHIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b8a3599d3f069ef59e0f0826b09e847d08d1d30d8dbe463188512a25a5a98be/data'
        : 'community.wave.seqera.io/library/cd-hit:4.8.1--2cbab405453b19db'}"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.fasta")    ,emit: fasta
    tuple val(meta), path("*.clstr")    ,emit: clusters
    tuple val("${task.process}"), val('cdhit'), eval("cd-hit -h | sed -n '1s/.*version \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_cdhit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[cd-hit] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    cd-hit \\
        $args \\
        -i ${sequences} \\
        -o ${prefix}.fasta \\
        -M $avail_mem \\
        -T $task.cpus
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.fasta
    touch ${prefix}.fasta.clstr
    """
}
