process HTSNIMTOOLS_VCFCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hts-nim-tools:0.3.11--hbeb723e_0':
        'biocontainers/hts-nim-tools:0.3.11--hbeb723e_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(background_vcf), path(background_tbi)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hts_nim_tools \\
        vcf-check \\
        $args \\
        $background_vcf \\
        $vcf \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htsnimtools: \$(hts_nim_tools | grep "version" | sed -e 's/version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htsnimtools: \$(hts_nim_tools | grep "version" | sed -e 's/version: //')
    END_VERSIONS
    """
}
