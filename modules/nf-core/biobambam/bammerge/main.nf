// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BIOBAMBAM_BAMMERGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4d04c55ec247a77fa53dba65675e35e714ef614cc549f2193058e14312d2448c/data'
        : 'community.wave.seqera.io/library/biobambam:2.0.185--f2942b2c30c9eb0a' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("*.bai"), optional: true, emit: bam_index
    tuple val(meta), path("*.md5"), optional: true, emit: checksum
    tuple val("${task.process}"), val('biobambam'), eval("bammerge --version |& sed '1!d; s/.*version //; s/.\$//'"), topic: versions, emit: versions_biobambam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_string = bam.join(" I=")

    """
    bammerge \\
        I=${input_string} \\
        ${args} \\
        > ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
