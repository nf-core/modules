process BISCUIT_BSCONV {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113" : null)
    def container_image = "biscuit:1.0.2.20220113--h81a5ba2_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(bam), path(bai)
    path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bsconv_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit bsconv \\
        $args \\
        \$INDEX \\
        $bam \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}
