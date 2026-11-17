// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/20/2009f48de55b9a633b784ce7acd3bd5d5d56fa74a0dc20543d8af624b6611faa/data'
        : 'community.wave.seqera.io/library/cancerit-allelecount:4.3.0--91e19f305ccc2f90'}"

    input:
    tuple val(meta), path(input), path(input_index)
    path loci
    path fasta

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    tuple val("${task.process}"), val('alleleCounter'), eval('alleleCounter --version'), emit: versions_allelecounter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_options = fasta ? "-r $fasta": ""
    """
    alleleCounter \\
        ${args} \\
        -l ${loci} \\
        -b ${input} \\
        ${reference_options} \\
        -o ${prefix}.alleleCount
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alleleCount
    """
}
