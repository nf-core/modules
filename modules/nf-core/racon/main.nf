// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RACON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6d/6da4caadf7fa541480c2498a763e7689350a97ee27ed8a661072ebcfd01f4bf9/data'
        : 'community.wave.seqera.io/library/racon:1.4.20--5a7230c7590a8e2e'}"

    input:
    tuple val(meta), path(reads), path(assembly), path(paf)

    output:
    tuple val(meta), path('*_assembly_consensus.fasta.gz') , emit: improved_assembly
    tuple val("${task.process}"), val('racon'), eval('racon --version 2>&1 | sed "s/^.*v//"'), emit: versions_racon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        "${paf}" \\
        $args \\
        "${assembly}" > \\
        ${prefix}_assembly_consensus.fasta

    gzip -n ${prefix}_assembly_consensus.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_assembly_consensus.fasta.gz
    """
}
