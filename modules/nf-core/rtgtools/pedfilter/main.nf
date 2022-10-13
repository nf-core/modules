process RTGTOOLS_PEDFILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::rtg-tools=3.12.1" : null)
    def container_image = "/rtg-tools:3.12.1--hdfd78af_0"
                                                     container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(ped)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rtg pedfilter \\
        $ped \\
        --vcf \\
        > ${prefix}.vcf

    gzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtgtools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
