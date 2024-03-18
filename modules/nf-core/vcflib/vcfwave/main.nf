process VCFLIB_VCFWAVE {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // WARN: vcfwave doesn't work in a conda environment (executable not found)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcflib:1.0.9--h146fbdb_4':
        'biocontainers/vcflib:1.0.9--h146fbdb_4' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    vcfwave \\
        $args \\
        --threads $task.cpus \\
        $vcf \\
        | bgzip $args --threads $task.cpus -c > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
