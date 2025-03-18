process MERQURYFK_KATGC {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/8948281b792d1d6385f4cf1f5da6d6139f04467ee73c171c4d54f43b23f4f4f4/data' :
        'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.2' }"

    input:
    tuple val(meta), path(fastk_hist), path(fastk_ktab)

    output:
    tuple val(meta), path("*.fi.{png,pdf}"), emit: filled_gc_plot , optional: true
    tuple val(meta), path("*.ln.{png,pdf}"), emit: line_gc_plot   , optional: true
    tuple val(meta), path("*.st.{png,pdf}"), emit: stacked_gc_plot, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    KatGC \\
        $args \\
        -T$task.cpus \\
        ${fastk_ktab.find{ it.toString().endsWith(".ktab") }} \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
