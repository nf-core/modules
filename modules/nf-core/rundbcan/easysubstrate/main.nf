// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RUNDBCAN_EASYSUBSTRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e843a691d1c9be41394116f30ee3cf3f7bd042e0d59cb2f6a0a34b2b0a308e30/data'
        : 'community.wave.seqera.io/library/dbcan:5.2.9--fb040f67fa175a9b'}"

    input:
    tuple val(meta), path(input_raw_data)
    tuple val(meta2), path(input_gff), val(gff_type)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}_overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_diamond.out"), emit: dbcandiamond_results
    tuple val(meta), path("${prefix}_cgc.gff"), emit: cgc_gff
    tuple val(meta), path("${prefix}_cgc_standard_out.tsv"), emit: cgc_standard_out
    tuple val(meta), path("${prefix}_diamond.out.tc"), emit: diamond_out_tc
    tuple val(meta), path("${prefix}_diamond.out.tf"), emit: diamond_out_tf, optional: true
    tuple val(meta), path("${prefix}_diamond.out.peptidase"), emit: diamond_out_peptidase
    tuple val(meta), path("${prefix}_diamond.out.sulfatase"), emit: diamond_out_sulfatase
    tuple val(meta), path("${prefix}_TF_hmm_results.tsv"), emit: tf_hmm_results, optional: true
    tuple val(meta), path("${prefix}_STP_hmm_results.tsv"), emit: stp_hmm_results
    tuple val(meta), path("${prefix}_total_cgc_info.tsv"), emit: total_cgc_info
    tuple val(meta), path("${prefix}_CGC.faa"), emit: cgc_faa
    tuple val(meta), path("${prefix}_PUL_blast.out"), emit: pul_blast
    tuple val(meta), path("${prefix}_substrate_prediction.tsv"), emit: substrate_prediction
    tuple val(meta), path("${prefix}_synteny_pdf/"), optional: true, emit: synteny_pdf
    tuple val("${task.process}"), val('rundbcan'), eval("run_dbcan version | sed 's/dbCAN version: //g'"), emit: versions_rundbcan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # substrate_prediction has no --threads and falls back to os.cpu_count()
    export PYTHON_CPU_COUNT=${task.cpus}

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir . \\
        --threads ${task.cpus} \\
        --input_gff ${input_gff} \\
        --gff_type ${gff_type} \\
        ${args}

    mv overview.tsv             ${prefix}_overview.tsv
    mv dbCAN_hmm_results.tsv    ${prefix}_dbCAN_hmm_results.tsv
    mv dbCANsub_hmm_results.tsv ${prefix}_dbCANsub_hmm_results.tsv
    mv diamond.out              ${prefix}_diamond.out
    mv cgc.gff                  ${prefix}_cgc.gff
    mv cgc_standard_out.tsv     ${prefix}_cgc_standard_out.tsv
    mv diamond.out.tc           ${prefix}_diamond.out.tc
    mv diamond.out.peptidase    ${prefix}_diamond.out.peptidase
    mv diamond.out.sulfatase    ${prefix}_diamond.out.sulfatase
    mv STP_hmm_results.tsv      ${prefix}_STP_hmm_results.tsv
    mv total_cgc_info.tsv       ${prefix}_total_cgc_info.tsv
    mv CGC.faa                  ${prefix}_CGC.faa
    mv PUL_blast.out            ${prefix}_PUL_blast.out
    mv substrate_prediction.tsv ${prefix}_substrate_prediction.tsv
    mv synteny_pdf/             ${prefix}_synteny_pdf/
    if [ -f diamond.out.tf ]; then
        mv diamond.out.tf       ${prefix}_diamond.out.tf
    fi
    if [ -f TF_hmm_results.tsv ]; then
        mv TF_hmm_results.tsv   ${prefix}_TF_hmm_results.tsv
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_tf_diamond = args.contains('--no-prokaryotic') ? '' : "touch ${prefix}_diamond.out.tf"
    def touch_tf_hmm     = args.contains('--fungi') ? "touch ${prefix}_TF_hmm_results.tsv" : ''
    """
    touch ${prefix}_overview.tsv
    touch ${prefix}_dbCAN_hmm_results.tsv
    touch ${prefix}_dbCANsub_hmm_results.tsv
    touch ${prefix}_diamond.out
    touch ${prefix}_cgc.gff
    touch ${prefix}_cgc_standard_out.tsv
    touch ${prefix}_diamond.out.tc
    ${touch_tf_diamond}
    touch ${prefix}_diamond.out.peptidase
    touch ${prefix}_diamond.out.sulfatase
    ${touch_tf_hmm}
    touch ${prefix}_STP_hmm_results.tsv
    touch ${prefix}_total_cgc_info.tsv
    touch ${prefix}_CGC.faa
    touch ${prefix}_PUL_blast.out
    touch ${prefix}_substrate_prediction.tsv
    mkdir -p ${prefix}_synteny_pdf
    """
}
