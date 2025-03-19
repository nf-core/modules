
process MERQURYFK_HAPMAKER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0bee9bfacd05665a9b1a11dd087dbf1be41ac3e640931c38c914a2390642cf/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-cowplot_r-ggplot2_r-viridis:f9994edc2270683c' }"

    input:
    tuple val(meta) , path(matktab)
    tuple val(meta2), path(patktab)
    tuple val(meta3), path(childktab)

    output:
    tuple val(meta) , path("*${input_mat}.hap.ktab*", hidden: true), emit: mathap_ktab
    tuple val(meta2), path("*${input_pat}.hap.ktab*", hidden: true), emit: pathap_ktab
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    input_mat = matktab ? "${matktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    input_pat = patktab ? "${patktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    def input_child = childktab ? "${childktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    HAPmaker \\
        $args \\
        -T${task.cpus} \\
        ${input_mat} \\
        ${input_pat} \\
        ${input_child}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FastK: $FASTK_VERSION
        MerquryFK: $MERQURY_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    input_mat = matktab ? "${matktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    input_pat = patktab ? "${patktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    touch ${input_mat}.hap.ktab
    touch .${input_mat}.hap.ktab.1
    touch ${input_pat}.hap.ktab
    touch .${input_pat}.hap.ktab.1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FastK: $FASTK_VERSION
        MerquryFK: $MERQURY_VERSION
    END_VERSIONS
    """
}