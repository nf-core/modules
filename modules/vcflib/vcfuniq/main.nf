def VERSION = '1.0.2' // Version information not provided by tool on CLI

process VCFLIB_VCFUNIQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcflib=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcflib:1.0.2--h3198e80_5':
        'quay.io/biocontainers/vcflib:1.0.2--h3198e80_5' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcfuniq \\
        $vcf \\
        | bgzip -c $args > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
