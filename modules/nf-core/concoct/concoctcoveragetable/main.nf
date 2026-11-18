// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CONCOCT_CONCOCTCOVERAGETABLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b414706a0ce09cd26055f920ea4d7db8e6e0a4845f806f0918e9792a2c847f95/data'
        : 'community.wave.seqera.io/library/concoct_setuptools:7d62fa8e6704b1fe'}"

    input:
    tuple val(meta), path(bed), path(bamfiles), path(baifiles)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('concoct'), eval('concoct --version | cut -d " " -f2'), emit: versions_concoct, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    concoct_coverage_table.py \\
        ${args} \\
        ${bed} \\
        ${bamfiles} \\
        > ${prefix}.tsv

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    """
}
