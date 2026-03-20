process SALTSHAKER_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c955cc086622ef50876a10e58a1e6711e42b70a0e4cbbc377142b62b0ad4f47/data':
        'community.wave.seqera.io/library/pip_saltshaker:ef543ea5ca09afbe' }"

    input:
    tuple val(meta), path(call)
    val dominant_fraction
    val group_radius
    val high_heteroplasmy
    val multiple_threshold
    val noise_threshold
    val mito_name

    output:
    tuple val(meta), path("*_classify_metadata.tsv"), emit: classify
    tuple val(meta), path("*_classify.txt")         , emit: txt
    tuple val(meta), path("*saltshaker.vcf")        , emit: vcf, optional: true
    tuple val("${task.process}"), val('saltshaker'), val("1.0.1"), topic: versions, emit: versions_saltshaker

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
        --chr-format $mito_name \\
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
