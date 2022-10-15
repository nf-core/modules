process VCFLIB_VCFBREAKMULTI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcflib=1.0.3" : null)
    def container_image = "vcflib:1.0.3--hecb563c_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vcfbreakmulti \\
        $vcf \\
        $args \\
        | bgzip -c $args2 > ${prefix}.breakmulti.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.breakmulti.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
