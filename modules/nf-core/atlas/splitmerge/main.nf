// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATLAS_SPLITMERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a494a4d166b29fdf9a27ac525c7aa4dc8ec476d0e624d174e6c72a0e87499dab/data'
        : 'community.wave.seqera.io/library/atlas:0.9.9--18077f78d3bd1a61'}"

    input:
    tuple val(meta), path(bam), path(bai), path(read_group_settings), path(blacklist)

    output:
    tuple val(meta), path("*_mergedReads.bam"), emit: bam
    tuple val(meta), path("*.txt.gz")         , emit: txt
    tuple val("${task.process}"), val('atlas'), eval("atlas | sed -e '2!d;s/.*Atlas //'"), emit: versions_atlas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def optional = blacklist ? 'blacklist=${blacklist}' : ''
    def args     = task.ext.args   ?: ''
    """
    atlas \\
        task=splitMerge bam=${bam} \\
        readGroupSettings=${read_group_settings}\\
        ${optional} \\
        ${args}
    """
}
