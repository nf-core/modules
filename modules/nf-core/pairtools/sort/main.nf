// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PAIRTOOLS_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc26df7ecf3557f443f89f4dc83380c68b413809cb6ca2213f14ea405c726080/data'
:         'community.wave.seqera.io/library/pairtools:1.1.3--f7f48429aa47aeaf' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: sorted
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*pairtools.*version //'") , emit: versions_pairtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer = task.memory.toGiga().intdiv(2)
    """
    pairtools \\
        sort \\
        $args \\
        --nproc $task.cpus \\
        --memory ${buffer}G \\
        -o ${prefix}.pairs.gz \\
        $input
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    """
}
