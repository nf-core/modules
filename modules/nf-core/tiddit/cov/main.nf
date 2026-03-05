process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89080e08e55823dbdb424407a4e9eefbc669e2b0e841f142a1014204659df87b/data' :
        'community.wave.seqera.io/library/tiddit:3.9.4--11a71ebb5fd55c20' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.bed"), optional: true, emit: cov
    tuple val(meta), path("${prefix}.wig"), optional: true, emit: wig
    tuple val("${task.process}"), val('tiddit'), eval("tiddit | sed -n 's/^usage: tiddit-//; s/ .*//p'"), topic: versions, emit: versions_tiddit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref ${fasta}" : ""
    """
    tiddit \\
        --cov \\
        -o ${prefix} \\
        ${args} \\
        --bam ${input} \\
        ${reference}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" > ${prefix}.wig
    echo "" > ${prefix}.bed
    """
}
