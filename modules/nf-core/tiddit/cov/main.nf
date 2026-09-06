process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3ebf66353bdc536851f786d7427e745c0c02f3951e08fc7828d6b86fefda64ce/data' :
        'community.wave.seqera.io/library/tiddit:3.9.7--33e2b6c6e2f37861' }"

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
