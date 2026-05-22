process FLUMUT_RUN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/flumut:0.6.5--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*xlsm"), emit: xlsm

    tuple val("${task.process}"), val('flumut'), eval("flumut --all-versions"), topic: versions, emit: versions_flumut

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flumut \\
        ${args} \\
        -x ${prefix}.xlsm \\
        -m ${prefix}_markers.tsv \\
        -M ${prefix}_mutations.tsv \\
        -l ${prefix}_literature.tsv \\
        ${fasta}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch \\
        ${prefix}.xlsm \\
        ${prefix}_markers.tsv \\
        ${prefix}_mutations.tsv \\
        ${prefix}_literature.tsv
    """
}
