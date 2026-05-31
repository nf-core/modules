process OPT_TRACK {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/nf-core/opt:v0.0.1"

    input:
    tuple val(meta), path(fwd_oriented_fa)
    tuple val(meta2), path(ref_annot_gff), path(ref_annot_fa)

    output:
    tuple val(meta), path("${prefix}/probe2targets.tsv"), emit: probes2target
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
        track \\
        -q ${fwd_oriented_fa} \\
        -a ${ref_annot_gff} \\
        -t ${ref_annot_fa} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/probe2targets.tsv"
    """
}
