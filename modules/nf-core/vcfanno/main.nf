// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VCFANNO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/ae01aba9334cde7a4a6f718930a7d158108332228f593140b52a2ec2727c336c/data'
        : 'community.wave.seqera.io/library/htslib_vcfanno:7ac3a954ec1a227e'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(specific_resources)
    path toml
    path lua
    path resources

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: tbi
    tuple val("${task.process}"), val('vcfanno'), eval("vcfanno 2>&1 | sed -n 's/.*version \\([0-9.]\\+\\).*/\\1/p'"), topic: versions, emit: versions_vcfanno

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def lua_cmd = lua ? "--lua ${lua}" : ""
    """
    vcfanno \\
        -p ${task.cpus} \\
        ${args} \\
        ${lua_cmd} \\
        ${toml} \\
        ${vcf} \\
        | bgzip ${args2} --threads ${task.cpus} \\
        > ${prefix}.vcf.gz

    tabix ${args3} ${prefix}.vcf.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    """
}
