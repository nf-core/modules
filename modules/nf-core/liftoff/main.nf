process LIFTOFF {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/liftoff:1.6.3--pyhdfd78af_0':
        'biocontainers/liftoff:1.6.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(target_fa)
    path ref_fa
    path ref_annotation
    path ref_db

    output:
    tuple val(meta), path("${prefix}.gff3")     , emit: gff3
    tuple val(meta), path("*.polished.gff3")    , emit: polished_gff3, optional: true
    tuple val(meta), path("*.unmapped.txt")     , emit: unmapped_txt
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?:  ''
    def arg_g   = ref_annotation    ?   "-g $ref_annotation"    : ''
    def arg_db  = ref_db            ?   "-db $ref_db"           : ''
    prefix      = task.ext.prefix   ?:  "${meta.id}"
    """
    if [[ ${target_fa} == *.gz ]]; then
        zcat ${target_fa} > target.fasta
    fi

    if [[ ${target_fa} == *.fa || ${target_fa} == *.fasta ]]; then
        ln -s ${target_fa} target.fasta
    fi

    if [[ ${ref_fa} == *.gz ]]; then
        zcat ${ref_fa} > reference.fasta
    fi

    if [[ ${ref_fa} == *.fa || ${ref_fa} == *.fasta ]]; then
        ln -s ${ref_fa} reference.fasta
    fi

    liftoff \\
        $arg_g \\
        $arg_db \\
        -p $task.cpus \\
        -o "${prefix}.gff3" \\
        -u "${prefix}.unmapped.txt" \\
        $args \\
        target.fasta \\
        reference.fasta

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
