// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b6f0f4d9e822f6def26b6f69a27cddd59475192c57c4a51385ef9ed13dd8078/data'
        : 'community.wave.seqera.io/library/methyldackel:0.6.1--aaf947cb078b2f30'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bedGraph") , optional: true, emit: bedgraph
    tuple val(meta), path("*.methylKit"), optional: true, emit: methylkit
    tuple val("${task.process}"), val('methyldackel'), eval("MethylDackel --version 2>&1 | cut -f1 -d' '"), emit: versions_methyldackel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    MethylDackel extract \\
        ${args} \\
        ${fasta} \\
        ${bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def out_extension = args.contains('--methylKit') ? 'methylKit' : 'bedGraph'
    """
    touch ${bam.baseName}_CpG.${out_extension}
    """
}
