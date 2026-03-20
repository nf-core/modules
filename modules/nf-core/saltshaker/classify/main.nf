process SALTSHAKER_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e93d703b195dd27cd920cee46669d3f51043216c12fd05168c937e93adf170e8/data':
        'community.wave.seqera.io/library/pip_saltshaker:e08e38a6d45f8f32' }"

    input:
    tuple val(meta), path(call)
    val dominant_fraction
    val group_radius
    val high_heteroplasmy
    val multiple_threshold
    val noise_threshold

    output:
    tuple val(meta), path("*_classify_metadata.tsv"), emit: classify
    tuple val(meta), path("*_classify.txt")         , emit: txt
    tuple val(meta), path("*saltshaker.vcf")        , emit: vcf, optional: true
    tuple val("${task.process}"), val('saltshaker'), val("1.0.0"), topic: versions, emit: versions_saltshaker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    saltshaker classify \\
        --prefix $prefix \\
        --input-dir . \\
        --dominant-fraction $dominant_fraction \\
        --radius $group_radius \\
        --high-het $high_heteroplasmy \\
        --multiple-threshold $multiple_threshold \\
        --noise $noise_threshold \\
        $args

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_vcf = args.contains('--vcf') ? "touch ${prefix}.saltshaker.vcf" : ''

    """
    echo $args

    $touch_vcf
    touch ${prefix}.saltshaker_classify.txt
    touch ${prefix}.saltshaker_classify_metadata.tsv
    """
}
