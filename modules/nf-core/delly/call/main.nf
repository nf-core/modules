process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::delly=1.1.5" : null)
        def container_image = "/delly:1.1.5--h358d541_0"
                                             container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(input), path(input_index), path(exclude_bed)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val(meta), path("*.csi"), emit: csi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    """
    delly \\
        call \\
        ${args} \\
        --outfile ${prefix}.bcf \\
        --genome ${fasta} \\
        ${exclude} \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
