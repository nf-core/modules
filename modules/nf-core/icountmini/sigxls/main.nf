process ICOUNTMINI_SIGXLS {
    tag "$meta.id"
    label "process_low"

    conda "bioconda::icount-mini=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0' :
        'biocontainers/icount-mini:2.0.3--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(bed)
    path segmentation

    output:
    tuple val(meta), path("*.sigxls.bed.gz"), emit: sigxls
    tuple val(meta), path("*.scores.tsv")   , emit: scores
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iCount-Mini sigxls \\
        $segmentation \\
        $bed \\
        ${prefix}.sigxls.bed.gz \\
        --scores ${prefix}.scores.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sigxls.bed.gz
    touch ${prefix}.scores.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}

