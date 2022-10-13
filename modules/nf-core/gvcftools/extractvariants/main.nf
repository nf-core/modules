process GVCFTOOLS_EXTRACTVARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gvcftools=0.17.0" : null)
        def container_image = "/gvcftools:0.17.0--he941832_3"
                                                            container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def open_gvcf = gvcf.extension == "gz" ? "gzip -dc $gvcf" : "cat $gvcf"

    """
    $open_gvcf |
    extract_variants \\
        $args \\
        $gvcf |
    gzip -c > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gvcftools: \$(extract_variants --help 2>&1 | grep version | sed 's/version: //')
    END_VERSIONS
    """
}
