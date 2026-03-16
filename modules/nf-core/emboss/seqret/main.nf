process EMBOSS_SEQRET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'biocontainers/emboss:6.6.0--h440b012_4' }"

    input:
    tuple val(meta), path(sequence)
    val out_ext

    output:
    tuple val(meta), path("*.${out_ext}"), emit: outseq
    tuple val("${task.process}"), val("emboss"), eval("revseq -version 2>&1 | sed 's/EMBOSS://'"), topic: versions, emit: versions_emboss

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.${out_ext}
    """
}
