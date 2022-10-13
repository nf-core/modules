process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    def container_image = "/lofreq:2.1.5--py38h588ecb2_4"
                                              container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq \\
        call \\
        $args \\
        -f $fasta \\
        -o ${prefix}.vcf.gz \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
