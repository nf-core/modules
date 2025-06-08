process CUSTOM_ADDMOSTSEVERECONSEQUENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbeb20c898a76bec809629320ece9e1f84a3e355e96568bfbe14b9f411bdf3e7/data':
        'community.wave.seqera.io/library/htslib_python:9c6265e98ef06930' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(variant_consequences)


    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: Version information not provided by tool on CLI.

    """
    add_most_severe_consequence.py \\
        $args \\
        --file_in $vcf \\
        --file_out ${prefix}.vcf \\
        --variant_csq $variant_consequences

    bgzip \\
        $args2 \\
        --threads ${task.cpus-1} \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence: $VERSION
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: Version information not provided by tool on CLI.
    """
    echo | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence: $VERSION
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """
}
