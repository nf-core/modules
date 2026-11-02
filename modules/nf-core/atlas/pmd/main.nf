// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ATLAS_PMD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a494a4d166b29fdf9a27ac525c7aa4dc8ec476d0e624d174e6c72a0e87499dab/data'
        : 'community.wave.seqera.io/library/atlas:0.9.9--18077f78d3bd1a61'}"

    input:
    tuple val(meta), path(bam), path(bai), path(pool_rg_txt)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*_PMD_input_Empiric.txt")    , emit: empiric
    tuple val(meta), path("*_PMD_input_Exponential.txt"), emit: exponential
    tuple val(meta), path("*_PMD_Table_counts.txt")     , emit: counts
    tuple val(meta), path("*_PMD_Table.txt")            , emit: table
    tuple val("${task.process}"), val('atlas'), eval('atlas | sed -e "2!d;s/.*Atlas //"'), emit: versions_atlas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def pool_rg_txt_cmd = pool_rg_txt ? "poolReadGroups=${pool_rg_txt}" : ""
    """
    atlas \\
        ${pool_rg_txt_cmd} \\
        task=PMD \\
        bam=${bam} \\
        fasta=${fasta} \\
        ${args}
    """
}
