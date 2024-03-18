process LIFTOFF {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/liftoff:1.6.3--pyhdfd78af_0':
        'biocontainers/liftoff:1.6.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(target_fa)
    path ref_fa, name: 'ref_assembly.fa'
    path ref_annotation

    output:
    tuple val(meta), path("${prefix}.gff3")     , emit: gff3
    tuple val(meta), path("*.polished.gff3")    , emit: polished_gff3, optional: true
    tuple val(meta), path("*.unmapped.txt")     , emit: unmapped_txt
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    liftoff \\
        -g $ref_annotation \\
        -p $task.cpus \\
        -o "${prefix}.gff3" \\
        -u "${prefix}.unmapped.txt" \\
        $args \\
        $target_fa \\
        ref_assembly.fa

    mv \\
        "${prefix}.gff3_polished" \\
        "${prefix}.polished.gff3" \\
        || echo "-polish is absent"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version 2> /dev/null)
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args     ?: ''
    prefix              = task.ext.prefix   ?: "${meta.id}"
    def touch_polished  = args.contains('-polish') ? "touch ${prefix}.polished.gff3" : ''
    """
    touch "${prefix}.gff3"
    touch "${prefix}.unmapped.txt"
    $touch_polished

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version 2> /dev/null)
    END_VERSIONS
    """
}
