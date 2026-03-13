
process SOMALIER_ANCESTRY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/somalier:0.2.19--h0c29559_0':
        'biocontainers/somalier:0.2.19--h0c29559_0' }"

    input:
    tuple val(meta),  path(query_somalier_files, stageAs: "query_files/*")
    tuple val(meta2), path(labels_tsv), path(labelled_somalier_files, stageAs: "labelled_files/*")

    output:
    tuple val(meta), path("*.somalier-ancestry.tsv") , emit: tsv
    tuple val(meta), path("*.somalier-ancestry.html"), emit: html
    tuple val("${task.process}"), val('somalier'), eval('somalier 2>&1 | sed -n \'s/.*version: \\([0-9.]*\\).*/\\1/p\''), emit: versions_somalier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    somalier  \\
        ancestry  \\
        --labels $labels_tsv \$(find -L labelled_files -type f -name "*.somalier")  \\
        ++ $query_somalier_files  \\
        --output-prefix=$prefix  \\
        $args
    """

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.somalier-ancestry.tsv"
    touch "${prefix}.somalier-ancestry.html"
    """

}
