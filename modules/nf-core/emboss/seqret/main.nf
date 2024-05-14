process EMBOSS_SEQRET {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::emboss=6.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'biocontainers/emboss:6.6.0--h440b012_4' }"

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
