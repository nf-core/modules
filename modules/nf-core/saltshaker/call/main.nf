process SALTSHAKER_CALL {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5a902cc9f161d602fde9c268a509be2f593cfac7ed4cdc2219f630e02e43b2ec/data':
        'community.wave.seqera.io/library/pip_saltshaker:be40ca61bbf77cf2' }"

    input:
    tuple val(meta), path(breakpoint), path(cluster)
    tuple val(meta2), path(mtfasta)
    val flank
    val heteroplasmy_limit
    val mito_length
    val heavy_strand_origin_start
    val heavy_strand_origin_end
    val light_strand_origin_start
    val light_strand_origin_end

    output:
    tuple val(meta), path("*_call_metadata.tsv"), emit: call
    tuple val("${task.process}"), val('saltshaker'), val("1.1.1"), topic: versions, emit: versions_saltshaker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    saltshaker call \\
        --prefix $prefix \\
        --output-dir . \\
        --reference $mtfasta \\
        --cluster $cluster \\
        --breakpoint $breakpoint \\
        --flank-size $flank \\
        --het-limit $heteroplasmy_limit \\
        --genome-length $mito_length \\
        --ori-h-start $heavy_strand_origin_start \\
        --ori-h-end $heavy_strand_origin_end \\
        --ori-l-start $light_strand_origin_start \\
        --ori-l-end $light_strand_origin_end \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.saltshaker_call_metadata.tsv
    """
}
