// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CIRCEXPLORER2_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/94b03d0007dd8092dfe7f057c6bdb521f772da234455e82112c12af3215af8ff/data'
        : 'community.wave.seqera.io/library/circexplorer2:2.3.8--2716c75300788f5b'}"

    input:
    tuple val(meta), path(fusions)

    output:
    tuple val(meta), path("*.bed"), emit: junction
    tuple val("${task.process}"), val('circexplorer2'), eval("CIRCexplorer2 --version 2>&1; true"), topic: versions, emit: versions_circexplorer2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligner = "${fusions}".endsWith(".junction") ? "-t STAR" : "${fusions}".endsWith(".txt") ? "-t MapSplice" : "${fusions}".endsWith(".bam") ? "-t BWA" : "-t segemehl"
    if ("${fusions}" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    CIRCexplorer2 \\
        parse \\
        $aligner \\
        $fusions \\
        -b ${prefix}.bed \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
