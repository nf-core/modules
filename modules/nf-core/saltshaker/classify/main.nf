process SALTSHAKER_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e93d703b195dd27cd920cee46669d3f51043216c12fd05168c937e93adf170e8/data':
        'community.wave.seqera.io/library/pip_saltshaker:e08e38a6d45f8f32' }"

    input:
    tuple val(meta), path(call)
    val dom_frac
    val group_radius
    val high_het
    val mult_thresh
    val noise_thresh

    output:
    tuple val(meta), path("*_classify_metadata.tsv"), emit: classify
    tuple val(meta), path("*_classify.txt")         , emit: txt
    tuple val(meta), path("*saltshaker.vcf")        , emit: vcf
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
        --vcf \\
        --dominant-fraction $dom_frac \\
        --radius $group_radius \\
        --high-het $high_het \\
        --multiple-threshold $mult_thresh \\
        --noise $noise_thresh \\
        $args

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.saltshaker.vcf
    touch ${prefix}.saltshaker_classify.txt
    touch ${prefix}.saltshaker_classify_metadata.tsv
    """
}
