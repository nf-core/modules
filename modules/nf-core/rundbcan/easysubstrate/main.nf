process RUNDBCAN_EASYSUBSTRATE {
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
    tuple val(meta), path("${prefix}/${prefix}_overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}/${prefix}_dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}/${prefix}_dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}/${prefix}_diamond.out"), emit: dbcandiamond_results
    tuple val(meta), path("${prefix}/${prefix}_cgc.gff"), emit: cgc_gff
    tuple val(meta), path("${prefix}/${prefix}_cgc_standard_out.tsv"), emit: cgc_standard_out
    tuple val(meta), path("${prefix}/${prefix}_diamond.out.tc"), emit: diamond_out_tc
    tuple val(meta), path("${prefix}/${prefix}_TF_hmm_results.tsv"), emit: tf_hmm_results
    tuple val(meta), path("${prefix}/${prefix}_STP_hmm_results.tsv"), emit: stp_hmm_results
    tuple val(meta), path("${prefix}/${prefix}_total_cgc_info.tsv"), emit: total_cgc_info
    tuple val(meta), path("${prefix}/${prefix}_CGC.faa"), emit: cgc_faa
    tuple val(meta), path("${prefix}/${prefix}_PUL_blast.out"), emit: pul_blast_out
    tuple val(meta), path("${prefix}/${prefix}_substrate_prediction.tsv"), emit: substrate_prediction
    tuple val(meta), path("${prefix}/${prefix}_synteny_pdf/"), emit: synteny_pdf




    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_substrate"
    def VERSION = '5.0.4'

    """
    mkdir -p ${prefix}
    touch ${prefix}/PUL_blast.out

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix} \\
        --input_gff ${input_gff} \\
        --gff_type ${gff_type} \\
        ${args}

    mv ${prefix}/overview.tsv             ${prefix}/${prefix}_overview.tsv
    mv ${prefix}/dbCAN_hmm_results.tsv    ${prefix}/${prefix}_dbCAN_hmm_results.tsv
    mv ${prefix}/dbCANsub_hmm_results.tsv ${prefix}/${prefix}_dbCANsub_hmm_results.tsv
    mv ${prefix}/diamond.out              ${prefix}/${prefix}_diamond.out
    mv ${prefix}/cgc.gff                  ${prefix}/${prefix}_cgc.gff
    mv ${prefix}/cgc_standard_out.tsv     ${prefix}/${prefix}_cgc_standard_out.tsv
    mv ${prefix}/diamond.out.tc           ${prefix}/${prefix}_diamond.out.tc
    mv ${prefix}/TF_hmm_results.tsv       ${prefix}/${prefix}_TF_hmm_results.tsv
    mv ${prefix}/STP_hmm_results.tsv      ${prefix}/${prefix}_STP_hmm_results.tsv
    mv ${prefix}/total_cgc_info.tsv       ${prefix}/${prefix}_total_cgc_info.tsv
    mv ${prefix}/CGC.faa                  ${prefix}/${prefix}_CGC.faa
    mv ${prefix}/PUL_blast.out            ${prefix}/${prefix}_PUL_blast.out
    mv ${prefix}/substrate_prediction.tsv ${prefix}/${prefix}_substrate_prediction.tsv
    mv ${prefix}/synteny_pdf/             ${prefix}/${prefix}_synteny_pdf/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_dbcan_substrate"
    def VERSION = '5.0.4'
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}_overview.tsv
    touch ${prefix}/${prefix}_dbCAN_hmm_results.tsv
    touch ${prefix}/${prefix}_dbCANsub_hmm_results.tsv
    touch ${prefix}/${prefix}_diamond.out
    touch ${prefix}/${prefix}_cgc.gff
    touch ${prefix}/${prefix}_cgc_standard_out.tsv
    touch ${prefix}/${prefix}_diamond.out.tc
    touch ${prefix}/${prefix}_TF_hmm_results.tsv
    touch ${prefix}/${prefix}_STP_hmm_results.tsv
    touch ${prefix}/${prefix}_total_cgc_info.tsv
    touch ${prefix}/${prefix}_CGC.faa
    touch ${prefix}/${prefix}_PUL_blast.out
    touch ${prefix}/${prefix}_substrate_prediction.tsv
    mkdir -p ${prefix}/${prefix}_synteny_pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """
}
