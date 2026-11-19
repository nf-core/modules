// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SVDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/17e4d7858e7b0aff0ed694e8a5f6955e8ba489fc3a9f0e9130ff2bf66d9ef1f0/data'
        : 'community.wave.seqera.io/library/htslib_svdss:10dd740771019782'}"

    input:
    tuple val(meta), path(bam), path(bai), path(sfs)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    SVDSS call \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        --bam ${bam} \\
        --sfs ${sfs} \\
        ${args} \\
        | bgzip -@ ${task.cpus} -c > ${prefix}.vcf.gz
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo "" | gzip > ${prefix}.vcf.gz
    """
}
