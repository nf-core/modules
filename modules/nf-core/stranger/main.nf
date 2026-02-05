process STRANGER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90f7f35b37e9e31c8680c18ef59c76e74500f44664887b3c1131daf07ff3043f/data':
        'community.wave.seqera.io/library/htslib_pip_stranger:de13e0cf88d77b50' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(variant_catalog)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('stranger'), eval("stranger --version | sed 's/stranger, version //g'"), topic: versions, emit: versions_stranger
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_stranger"
    def options_variant_catalog = variant_catalog ? "--repeats-file $variant_catalog" : ""

    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    stranger \\
        $args \\
        $vcf \\
        $options_variant_catalog | bgzip $args2 -c --threads ${task.cpus} > ${prefix}.vcf.gz

    tabix \\
        $args3 \\
        --threads ${task.cpus} \\
        ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_stranger"

    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
