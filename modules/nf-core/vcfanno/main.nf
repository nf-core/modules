process VCFANNO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b23fffa38d9740616e3b414df24e44c22fbc0510264f79b5062b3eecd619393f/data':
        'community.wave.seqera.io/library/htslib_vcfanno:c88d6077509197fe' }"

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
        > ${prefix}.vcf.gz \\
        && tabix ${args3} ${prefix}.vcf.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    """
}
