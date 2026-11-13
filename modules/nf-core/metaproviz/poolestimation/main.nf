process METAPROVIZ_POOLESTIMATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-metaproviz:4.99.0--r45hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-metaproviz:4.99.0--r45hdfd78af_0' }"

    input:
    tuple val(meta), path(se_rds), path(data_matrix), path(feature_matrix), path(sample_matrix)
    val output_format

    output:
    tuple val(meta), path("*.pool_filtered.rds")                       , optional: true, emit: se
    tuple val(meta), path("*.pool_filtered.*_matrix.{tsv,csv}")        , optional: true, emit: filtered_matrices
    tuple val(meta), path("*.pools.txt")                               , emit: pools
    tuple val(meta), path("*.cv.tsv")                                  , emit: cv_table
    tuple val(meta), path("*.plots.rds")                               , emit: plots
    tuple val(meta), path("*.report.html")                             , emit: report
    tuple val(meta), path("*.log")                                     , emit: log
    tuple val("${task.process}"), val('metaproviz'), val('4.99.0'), emit: versions_metaproviz, topic: versions
    tuple val("${task.process}"), val('r-base'), val('4.5.3'), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'metaproviz_poolestimation.R'

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pool_filtered.rds ${prefix}.pools.txt ${prefix}.cv.tsv \\
          ${prefix}.pool_filtered.{data_matrix,feature_matrix,sample_matrix}.tsv \\
          ${prefix}.plots.rds ${prefix}.report.html ${prefix}.log
    """
}
