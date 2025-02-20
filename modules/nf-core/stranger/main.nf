process STRANGER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stranger:0.9.2--pyh7e72e81_0':
        'biocontainers/stranger:0.9.2--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(variant_catalog)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_variant_catalog = variant_catalog ? "--repeats-file $variant_catalog" : ""
    """
    stranger \\
        $args \\
        $vcf \\
        $options_variant_catalog | gzip --no-name > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$( stranger --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$( stranger --version )
    END_VERSIONS
    """
}
