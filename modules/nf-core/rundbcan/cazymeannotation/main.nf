process RUNDBCAN_CAZYMEANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.0.4--pyhdfd78af_0' :
        'biocontainers/dbcan:5.0.4--pyhdfd78af_0' }"
    input:
    tuple val(meta), path(input_raw_data)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}/overview.tsv")            , emit: cazyme_annotation
    tuple val(meta), path("${prefix}/dbCAN_hmm_results.tsv")   , emit: dbcanhmm_results
    tuple val(meta), path("${prefix}/dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}/diamond.out")             , emit: dbcandiamond_results

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_cazyme"
    def VERSION = '5.0.4'

    """

    run_dbcan CAZyme_annotation \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_cazyme"
    def VERSION = '5.0.4'
    """
    mkdir -p ${prefix}
    touch ${prefix}/overview.tsv
    touch ${prefix}/dbCAN_hmm_results.tsv
    touch ${prefix}/dbCANsub_hmm_results.tsv
    touch ${prefix}/diamond.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """
}
