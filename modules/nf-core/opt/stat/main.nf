process OPT_STAT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/nf-core/opt:v0.0.1"

    input:
    tuple val(meta), path(probe_targets)
    tuple val(meta2), path(fwd_oriented_probes)
    path(gene_synonyms)

    output:
    tuple val(meta), path("${prefix}/collapsed_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('opt'), eval("opt --version"), topic: versions, emit: versions_opt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def synonyms = gene_synonyms ? "-s ${gene_synonyms}": ""

    """
    opt \\
        -o ${prefix} \\
        stat \\
        -i ${probe_targets} \\
        -q ${fwd_oriented_probes} \\
        ${synonyms} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/collapsed_summary.tsv"
    """
}
