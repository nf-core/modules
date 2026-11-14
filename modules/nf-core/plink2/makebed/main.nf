// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLINK2_MAKEBED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d0745d58f97f7ada3bad4ece2d69f97922da71c47cb4530cce8d7264b8142990/data'
        : 'community.wave.seqera.io/library/plink2:2.0.0a7.post13--1df32678f52384f3'}"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)

    output:
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam"), emit: bed
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//; s/LM\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --pgen ${pgen} \\
        --pvar ${pvar} \\
        --psam ${psam} \\
        --threads ${task.cpus} \\
        --memory ${mem_mb} \\
        ${args} \\
        --make-bed \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bim
    touch ${prefix}.fam
    touch ${prefix}.log
    """
}
