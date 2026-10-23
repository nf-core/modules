// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BAMSTATS_GENERALSTATS {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/74dd96801e132db6cada326a68fc9c6e1e3fe08cf2e3846787d1d4ff29433982/data'
        : 'community.wave.seqera.io/library/bamstats:0.3.5--0a6b61f4dd590324' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val("bamstats"), eval('bamstats --version | grep "version: " | sed -e s"/version: //"'), topic: versions, emit: versions_bamstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Several ARGS are available.
    //  -a is a helpfu; one where you can add a BED file
    //  -u, --uniq outputs genomic coverage statistics for uniqely mapped reads
    """
    bamstats \\
        -i ${bam} \\
        ${args} \\
        -c ${task.cpus} \\
        -o ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.json
    """
}
