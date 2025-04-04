process RUNDBCAN_EASYSUBSTRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/bcb-unl/run_dbcan_new:5.0.2" // Use the same container as RUNDBCAN_DATABASE

    input:
    tuple val(meta), path(input_raw_data)
    tuple val(meta), path(input_gff)
    tuple val(meta), val (gff_type)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}_dbcan_substrate/overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbcan_substrate/dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbcan_substrate/dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_dbcan_substrate/diamond.out"), emit: dbcandiamond_results
    tuple val(meta), path("${prefix}_dbcan_substrate/cgc.gff"), emit: cgc_gff
    tuple val(meta), path("${prefix}_dbcan_substrate/cgc_standard_out.tsv"), emit: cgc_standard_out
    tuple val(meta), path("${prefix}_dbcan_substrate/diamond.out.tc"), emit: diamond_out_tc
    tuple val(meta), path("${prefix}_dbcan_substrate/TF_hmm_results.tsv"), emit: tf_hmm_results
    tuple val(meta), path("${prefix}_dbcan_substrate/STP_hmm_results.tsv"), emit: stp_hmm_results
    tuple val(meta), path("${prefix}_dbcan_substrate/total_cgc_info.tsv"), emit: total_cgc_info
    tuple val(meta), path("${prefix}_dbcan_substrate/CGC.faa"), emit: cgc_faa
    tuple val(meta), path("${prefix}_dbcan_substrate/PUL_blast.out"), emit: pul_blast_out
    tuple val(meta), path("${prefix}_dbcan_substrate/substrate_prediction.tsv"), emit: substrate_prediction
    tuple val(meta), path("${prefix}_dbcan_substrate/synteny_pdf/"), emit: synteny_pdf




    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.0.2'

    """

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix}_dbcan_substrate \\
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
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.0.2'
    """
    mkdir -p ${prefix}_dbcan_substrate
    touch ${prefix}_dbcan_substrate/overview.tsv
    touch ${prefix}_dbcan_substrate/dbCAN_hmm_results.tsv
    touch ${prefix}_dbcan_substrate/dbCANsub_hmm_results.tsv
    touch ${prefix}_dbcan_substrate/diamond.out
    touch ${prefix}_dbcan_substrate/cgc.gff
    touch ${prefix}_dbcan_substrate/cgc_standard_out.tsv
    touch ${prefix}_dbcan_substrate/diamond.out.tc
    touch ${prefix}_dbcan_substrate/TF_hmm_results.tsv
    touch ${prefix}_dbcan_substrate/STP_hmm_results.tsv
    touch ${prefix}_dbcan_substrate/total_cgc_info.tsv
    touch ${prefix}_dbcan_substrate/CGC.faa
    touch ${prefix}_dbcan_substrate/PUL_blast.out
    touch ${prefix}_dbcan_substrate/substrate_prediction.tsv
    mkdir -p ${prefix}_dbcan_substrate/synteny_pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """
}
