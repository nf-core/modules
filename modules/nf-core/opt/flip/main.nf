process OPT_FLIP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/khersameesh24/opt:v0.0.1"

    input:
    tuple val(meta), path(probes_fasta)
    tuple val(meta2), path(ref_annot_gff), path(ref_annot_fa)

    output:
    tuple val(meta), path("${prefix}/fwd_oriented.fa"), emit: fwd_oriented_fa
    tuple val("${task.process}"), val('opt'), eval("opt --version"), topic: versions, emit: versions_opt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    opt \\
        -o ${prefix} \\
        -p ${task.cpus} \\
        flip \\
        -q ${probes_fasta} \\
        -a ${ref_annot_gff} \\
        -t ${ref_annot_fa} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/fwd_oriented.fa"
    """
}
