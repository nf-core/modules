process LOFREQ_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    def container_image = "lofreq:2.1.5--py38h588ecb2_4"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq \\
        filter \\
        $args \\
        -i $vcf \\
        -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
