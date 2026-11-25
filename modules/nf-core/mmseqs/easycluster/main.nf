// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_EASYCLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/161e2ded0fc1b8d724333808695d7da12a8643668a59b652796e7a51ff9b4f3e/data'
        : 'community.wave.seqera.io/library/mmseqs2:20.ae60e--0746edfc81192f32'}"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*rep_seq.fasta"), emit: representatives
    tuple val(meta), path("*all_seqs.fasta"), emit: fasta
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mmseqs \\
        easy-cluster \\
        ${sequence} \\
        ${prefix} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.tsv
    touch ${prefix}_rep_seq.fasta
    touch ${prefix}_all_seqs.fasta

    """
}
