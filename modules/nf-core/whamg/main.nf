// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process WHAMG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b5cba56edb5c342cde8d5747311d90b4678391ee383bbfbf74dd2d524e2cbcc6/data'
:         'community.wave.seqera.io/library/tabix_wham:0eb67c2c9cb26529' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: tbi
    tuple val(meta), path("*.txt")          , emit: graph, optional: true
    tuple val("${task.process}"), val('whamg'), eval("whamg 2>&1 | sed -n 's/^Version: v\\([^-]*\\).*/\\1/p'"), topic: versions, emit: versions_whamg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    whamg \\
        -x ${task.cpus} \\
        -f ${bam} \\
        -a ${fasta} \\
        ${args} \\
        | bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz

        tabix ${args3} ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
