// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CLAME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b4459c82bcfba35249d0db055fe5b8d37c15a703052b52e9c2bfeb2b7c1a724d/data'
        : 'community.wave.seqera.io/library/clame:1.0--a3c510a4b70a1afd'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta")  , emit: fasta, optional: true
    tuple val(meta), path("*.binning"), emit: bins
    tuple val(meta), path("*.fm9")    , emit: fm
    tuple val(meta), path("*.index")  , emit: index
    tuple val(meta), path("*.links")  , emit: links
    tuple val(meta), path("*.result") , emit: result
    tuple val("${task.process}"), val('clame'), eval("clame -h 2>&1 | sed '2!d;s/version //;s/ .*//'"), emit: versions_clame, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clame \\
        $args \\
        -nt $task.cpus \\
        -multiFasta ${fasta} \\
        -output ${prefix} || test -f ${prefix}.binning
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.binning
    touch ${prefix}.fm9
    touch ${prefix}.index
    touch ${prefix}.links
    touch ${prefix}.result
    """
}
