process POLYPOLISH_POLISH {

    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/polypolish:0.6.0--hdbdd923_0':
        'biocontainers/polypolish:0.6.0--hdbdd923_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(sam)
    val save_debug

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.txt"), optional: true, emit: debug
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def debug_mode = save_debug ? "--debug ${prefix}.txt" : ''

    """
    polypolish \\
        polish \\
        $args \\
        $debug_mode \\
        $fasta \\
        $sam > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish polish --version |& sed '1!d ; s/Polypolish-polish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish polish --version |& sed '1!d ; s/polypolish polish //')
    END_VERSIONS
    """
}
