
process SOMALIER_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::somalier=0.2.15" : null)
        def container_image = "/somalier:0.2.15--h37c5b7d_0"
                                                   container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(ref), path(refidx)
    path(sites)

    output:
    tuple val(meta), path("*.somalier"), emit: extract
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    somalier extract \\
    --sites ${sites} \\
    -f ${ref} \\
    ${bam} \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """
}
