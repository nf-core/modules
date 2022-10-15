process EMBOSS_SEQRET {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)
    def container_image = "emboss:6.6.0--hf657eab_5"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(sequence)
    val out_ext

    output:
    tuple val(meta), path("*.${out_ext}"), emit: outseq
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def osformat = args.contains('-osformat') ? '' : "-osformat ${out_ext}"
    """
    seqret \\
        ${args} \\
        -sequence ${sequence} \\
        ${osformat} \\
        -outseq ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(echo \$(seqret -version 2>&1) | sed 's/EMBOSS://')
    END_VERSIONS
    """
}
