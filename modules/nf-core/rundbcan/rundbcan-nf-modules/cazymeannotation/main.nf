process RUNDBCAN_CAZYMEANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/bcb-unl/run_dbcan_new:5.0.2" // Use the same container as RUNDBCAN_DATABASE

    input:
    tuple val(meta), path(input_raw_data)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}_dbcan_cazyme/overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbcan_cazyme/dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbcan_cazyme/dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_dbcan_cazyme/diamond.out"), emit: dbcandiamond_results

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.0.2'

    """

    run_dbcan CAZyme_annotation \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix}_dbcan_cazyme \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.0.2'
    """
    mkdir -p ${prefix}_dbcan_cazyme
    touch ${prefix}_dbcan_cazyme/overview.tsv
    touch ${prefix}_dbcan_cazyme/dbCAN_hmm_results.tsv
    touch ${prefix}_dbcan_cazyme/dbCANsub_hmm_results.tsv
    touch ${prefix}_dbcan_cazyme/diamond.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """
}
