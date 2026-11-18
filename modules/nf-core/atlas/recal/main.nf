// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATLAS_RECAL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a494a4d166b29fdf9a27ac525c7aa4dc8ec476d0e624d174e6c72a0e87499dab/data'
        : 'community.wave.seqera.io/library/atlas:0.9.9--18077f78d3bd1a61'}"

    input:
    tuple val(meta), path(bam), path(bai), path(empiric), path(readgroups)
    path(alleles)
    path(invariant_sites)

    output:
    tuple val(meta), path("*.txt"), emit:recal_patterns
    tuple val("${task.process}"), val('atlas'), eval('((atlas 2>&1) | grep Atlas | head -n 1 | sed -e \'s/^[ \t]*Atlas //\')'), emit: versions_atlas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def PMD        = empiric         ? "pmdFile=${empiric}"           : ""
    def ALLELES    = alleles         ? "alleleFile=${alleles}"        : ""
    def INVARIANTS = invariant_sites ? "window=${invariant_sites}"    : ""
    def READGROUPS = readgroups      ? "poolReadGroups=${readgroups}" : ""

    """
    atlas \\
        task=recal \\
        bam=${bam} \\
        ${PMD} \\
        ${READGROUPS} \\
        ${ALLELES} \\
        ${INVARIANTS} \\
        out=${prefix} \\
        ${args}

    """
}
