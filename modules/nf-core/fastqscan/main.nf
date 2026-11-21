// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTQSCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e87ef091f75815b874e21b0923c667168cae5e3443b6d5a469dc30b572a669c/data'
        : 'community.wave.seqera.io/library/fastq-scan:1.0.1--107ab5ca81ecea03'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val('fastqscan'), eval('fastq-scan -v 2>&1 | sed \'s/^.*fastq-scan //\''), emit: versions_fastqscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $reads | \\
        fastq-scan \\
        $args > ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    """
}
