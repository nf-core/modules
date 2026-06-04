process ICOUNTMINI_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:3.0.1--pyh7cba7a3_0':
        'quay.io/biocontainers/icount-mini:3.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bed)
    path segmentation

    output:
    tuple val(meta), path("*summary_type.tsv")   , emit: summary_type
    tuple val(meta), path("*summary_subtype.tsv"), emit: summary_subtype
    tuple val(meta), path("*summary_gene.tsv")   , emit: summary_gene
    tuple val("${task.process}"), val('iCount-Mini'), eval("iCount-Mini -v"), emit: versions_icount_mini, topic: versions


    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini summary \\
        ${segmentation} \\
        ${bed} \\
        . \\
        ${args}

    mv summary_type.tsv ${prefix}.summary_type.tsv
    mv summary_subtype.tsv ${prefix}.summary_subtype.tsv
    mv summary_gene.tsv ${prefix}.summary_gene.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary_type.tsv
    touch ${prefix}.summary_subtype.tsv
    touch ${prefix}.summary_gene.tsv
    """
}
