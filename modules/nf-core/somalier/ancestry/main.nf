
process SOMALIER_ANCESTRY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::somalier=0.2.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/somalier:0.2.15--h37c5b7d_0':
        'quay.io/biocontainers/somalier:0.2.15--h37c5b7d_0' }"

    input:
    tuple val(meta), path(query_somalier_files)
    path(labels_tsv)
    path(labelled_somalier_files_directory)

    output:
    tuple val(meta), path("*-ancestry.tsv"), emit: tsv
    tuple val(meta), path("*-ancestry.html"), emit: html
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    labelled_somalier_files = "${labelled_somalier_files_directory}/*.somalier"

    """
    somalier  \\
        ancestry  \\
        --labels $labels_tsv $labelled_somalier_files  \\
        ++ $query_somalier_files  \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """

}
