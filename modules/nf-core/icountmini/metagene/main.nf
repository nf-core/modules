process ICOUNTMINI_METAGENE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::icount-mini=3.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:3.0.1--pyh7cba7a3_0':
        'biocontainers/icount-mini:3.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bed)
    path segmentation

    output:
    tuple val(meta), path("metagene_*/*plot_data.tsv"), emit: tsv
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mv $bed ${prefix}.bed

    iCount-Mini metagene \\
        ${prefix}.bed \\
        $segmentation \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch metagene_${prefix}/${prefix}_plot_data.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}
