// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process EPANG_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5227fa1db763e3ed8267ff4f9e48d02b6a8a595a065607fd1a8f0b7fb8b17bb1/data'
        : 'community.wave.seqera.io/library/epa-ng:0.3.8--b7414cbf4b37d776'}"

    input:
    tuple val(meta), path(refaln), path(fullaln)

    output:
    tuple val(meta), path("*query.fasta.gz")    , emit: query
    tuple val(meta), path("*reference.fasta.gz"), emit: reference
    tuple val("${task.process}"), val('epa-ng'), eval('epa-ng --version | sed "s/EPA-ng v//"'), emit: versions_epang, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epa-ng \\
        $args \\
        --split $refaln $fullaln

    gzip -c query.fasta > ${prefix}.query.fasta.gz; rm query.fasta
    gzip -c reference.fasta > ${prefix}.reference.fasta.gz; rm reference.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_query.fasta
    gzip ${prefix}_query.fasta
    touch ${prefix}_reference.fasta
    gzip ${prefix}_reference.fasta
    """
}
