process SUPPA_DIFFSPLICE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    tuple val(meta), path(events)
    tuple val(meta2), val(condition1), path(expression1), path(psi1)
    tuple val(meta3), val(condition2), path(expression2), path(psi2)
    val method
    val local_area
    val lower_bound_delta_psi
    val is_paired
    val gene_correction
    val multi_testing_correction_alpha
    val save_tpm_events
    val combination_analysis
    val use_median_delta_psi
    val tpm_threshold
    val nan_tpm_threshold

    output:
    tuple val(meta), path("*.dpsi"), emit: dpsi
    tuple val(meta), path("*.psivec"), emit: psivec
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_arg = is_paired ? '--paired' : ''
    def area_arg = local_area ? "--area ${local_area}" : ''
    def lower_bound_delta_psi_arg = lower_bound_delta_psi ? "--lower-bound ${lower_bound_delta_psi}" : ''
    def gene_correction_arg = gene_correction ? '--gene-correction' : ''
    def combination_analysis_arg = combination_analysis ? '--combination' : ''
    def alpha_arg = multi_testing_correction_alpha ? "--alpha ${multi_testing_correction_alpha}" : ''
    def median_delta_psi_arg = use_median_delta_psi ? '--median' : ''
    def tpm_threshold_arg = tpm_threshold ? "--tpm-threshold ${tpm_threshold}" : ''
    def nan_tpm_threshold_arg = nan_tpm_threshold ? "--nan-tpm-threshold ${nan_tpm_threshold}" : ''
    """
    suppa.py \\
        diffSplice \\
        --method ${method} \\
        --psi ${psi1} ${psi2} \\
        --tpm ${expression1} ${expression2} \\
        --input ${events} \\
        --output ${prefix} \\
        ${area_arg} \\
        ${lower_bound_delta_psi_arg} \\
        ${paired_arg} \\
        ${gene_correction_arg} \\
        ${combination_analysis_arg} \\
        ${alpha_arg} \\
        ${median_delta_psi_arg} \\
        ${tpm_threshold_arg} \\
        ${nan_tpm_threshold_arg} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.dpsi
    touch ${prefix}.psivec
    """
}
