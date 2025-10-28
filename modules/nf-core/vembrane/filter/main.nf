process VEMBRANE_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vembrane:2.4.0--pyhdfd78af_0':
        'biocontainers/vembrane:2.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(variant)
    val(expression)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: filtered_variant
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vembrane filter \\
        ${args} \\
        ${expression} \\
        -o ${prefix}_filtered.vcf \\
        $variant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vembrane: \$(vembrane --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    vembrane filter \\
        ${args} \\
        ${expression} \\
        -o ${prefix}_filtered.vcf \\
        $variant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vembrane: \$(vembrane --version)
    END_VERSIONS
    """
}
