// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEEPTOOLS_BIGWIGCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba9641d61acf134c6b989b3af5f25bcd5aadb3faf692a41c00ba4cd25f99032e/data'
:         'community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3' }"

    input:
    tuple val(meta) , path(bigwig1)     , path(bigwig2)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.{bigWig,bedgraph}"), emit: output
    tuple val("${task.process}"), val('deeptools'), eval('bigwigCompare --version | sed "s/bigwigCompare //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args                                  ?: ""
    def prefix = task.ext.prefix                              ?: "${meta.id}"
    def blacklist_cmd = blacklist                             ? "--blackListFileName ${blacklist}" : ""
    def extension = args.contains("--outFileFormat bedgraph") ? "bedgraph"                         : "bigWig"

    """
    bigwigCompare \\
        --bigwig1 $bigwig1 \\
        --bigwig2 $bigwig2 \\
        --outFileName ${prefix}.${extension} \\
        --numberOfProcessors $task.cpus \\
        $blacklist_cmd \\
        $args
    """

    stub:
    def args = task.ext.args                                  ?: ''
    def prefix = task.ext.prefix                              ?: "${meta.id}"
    def extension = args.contains("--outFileFormat bedgraph") ? "bedgraph" : "bigWig"

    """
    touch ${prefix}.${extension}
    """
}
