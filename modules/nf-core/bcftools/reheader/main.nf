process BCFTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    def container_image = "bcftools:1.15.1--h0ea216a_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(vcf)
    path fai
    path header

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def update_sequences = fai ? "-f $fai" : ""
    def new_header       = header ? "-h $header" : ""
    """
    bcftools \\
        reheader \\
        $update_sequences \\
        $new_header \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix}.vcf.gz \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
