process ICOUNTMINI_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::icount-mini=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0':
        'biocontainers/icount-mini:2.0.3--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(bed)
    path segmentation

    output:
    tuple val(meta), path("*summary_type.tsv")   , emit: summary_type
    tuple val(meta), path("*summary_subtype.tsv"), emit: summary_subtype
    tuple val(meta), path("*summary_gene.tsv")   , emit: summary_gene
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini summary \\
        $segmentation \\
        $bed \\
        . \\
        $args

    mv summary_type.tsv ${prefix}.summary_type.tsv
    mv summary_subtype.tsv ${prefix}.summary_subtype.tsv
    mv summary_gene.tsv ${prefix}.summary_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary_type.tsv
    touch ${prefix}.summary_subtype.tsv
    touch ${prefix}.summary_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}
