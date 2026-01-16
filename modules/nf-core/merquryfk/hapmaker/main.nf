
process MERQURYFK_HAPMAKER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e7bf3e855f8b4abb82024913359efef5dfc40902913bede1a6e6fbf8fc110bc/data' :
        'community.wave.seqera.io/library/fastk_merquryfk_r-ggplot2_r-viridis_pruned:bbe5238acf370901' }"

    input:
    tuple val(meta) , path(matktab)
    tuple val(meta2), path(patktab)
    tuple val(meta3), path(childktab)

    output:
    tuple val(meta) , path("*${input_mat}.hap.ktab*", hidden: true), emit: mat_hap_ktab
    tuple val(meta2), path("*${input_pat}.hap.ktab*", hidden: true), emit: pat_hap_ktab
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('merquryfk'), eval('echo 1.1.3'), emit: versions_merquryfk, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastk'), eval('echo 1.1'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    input_mat       = matktab   ? "${matktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}"   : ''
    input_pat       = patktab   ? "${patktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}"   : ''
    def input_child = childktab ? "${childktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    """
    HAPmaker \\
        $args \\
        -T${task.cpus} \\
        ${input_mat} \\
        ${input_pat} \\
        ${input_child}
    """

    stub:
    def args  = task.ext.args ?: ''
    input_mat = matktab ? "${matktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    input_pat = patktab ? "${patktab.find{ it.toString().endsWith(".ktab") }.toString() - ~/\.ktab/}" : ''
    """
    touch ${input_mat}.hap.ktab
    touch .${input_mat}.hap.ktab.1
    touch ${input_pat}.hap.ktab
    touch .${input_pat}.hap.ktab.1
    """
}
