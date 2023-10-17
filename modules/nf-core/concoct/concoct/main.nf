
process CONCOCT_CONCOCT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::concoct=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py311h245ed52_4':
        'biocontainers/concoct:1.1.0--py311h245ed52_4' }"

    input:
    tuple val(meta), path(coverage_file), path(fasta)

    output:
    tuple val(meta), path("*_args.txt")                         , emit: args_txt
    tuple val(meta), path("*_clustering_gt1000.csv")            , emit: clustering_csv
    tuple val(meta), path("*_log.txt")                          , emit: log_txt
    tuple val(meta), path("*_original_data_gt1000.csv")         , emit: original_data_csv
    tuple val(meta), path("*_PCA_components_data_gt1000.csv")   , emit: pca_components_csv
    tuple val(meta), path("*_PCA_transformed_data_gt1000.csv")  , emit: pca_transformed_csv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    concoct \\
        $args \\
        --threads ${task.cpus} \\
        --coverage_file ${coverage_file} \\
        --composition_file ${fasta} \\
        -b ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """
}
