process STRANGER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/19f14a7c1b0ec9cbdeb6d32e3692208d559e9186b210e9a0a6922e001cb6ad32/data':
        'community.wave.seqera.io/library/stranger_tabix:847b205e87ed124b' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(variant_catalog)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$( stranger --version | sed 's/stranger, version //g' )
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_stranger"

    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$( stranger --version | sed 's/stranger, version //g' )
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
