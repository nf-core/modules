process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0bee9bfacd05665a9b1a11dd087dbf1be41ac3e640931c38c914a2390642cf/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-cowplot_r-ggplot2_r-viridis:f9994edc2270683c' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*.hist"), emit: hist
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    Histex \\
        $args \\
        $histogram \\
        > ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
