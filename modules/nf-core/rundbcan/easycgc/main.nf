process RUNDBCAN_EASYCGC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.0.4--pyhdfd78af_0' :
        'biocontainers/dbcan:5.0.4--pyhdfd78af_0' }"
    input:
    tuple val(meta), path(input_raw_data)
    tuple val(meta2), path(input_gff), val (gff_type)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}/overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}/dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}/dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}/diamond.out"), emit: dbcandiamond_results
    tuple val(meta), path("${prefix}/cgc.gff"), emit: cgc_gff
    tuple val(meta), path("${prefix}/cgc_standard_out.tsv"), emit: cgc_standard_out
    tuple val(meta), path("${prefix}/diamond.out.tc"), emit: diamond_out_tc
    tuple val(meta), path("${prefix}/TF_hmm_results.tsv"), emit: tf_hmm_results
    tuple val(meta), path("${prefix}/STP_hmm_results.tsv"), emit: stp_hmm_results
    tuple val(meta), path("${prefix}/total_cgc_info.tsv"), emit: total_cgc_info

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_cgc"
    def VERSION = '5.0.4'

    """

    run_dbcan easy_CGC \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix} \\
        --input_gff ${input_gff} \\
        --gff_type ${gff_type} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_cgc"
    def VERSION = '5.0.4'
    """
    mkdir -p ${prefix}
    touch ${prefix}/overview.tsv
    touch ${prefix}/dbCAN_hmm_results.tsv
    touch ${prefix}/dbCANsub_hmm_results.tsv
    touch ${prefix}/diamond.out
    touch ${prefix}/cgc.gff
    touch ${prefix}/cgc_standard_out.tsv
    touch ${prefix}/diamond.out.tc
    touch ${prefix}/TF_hmm_results.tsv
    touch ${prefix}/STP_hmm_results.tsv
    touch ${prefix}/total_cgc_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """
}
