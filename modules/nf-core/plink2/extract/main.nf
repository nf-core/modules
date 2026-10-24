// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLINK2_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/278e4e95ad25ec7048ab124f16a7cf8b8141ced9e05c51f61803567cd9017f4d/data'
:         'community.wave.seqera.io/library/plink2:2.00a2.3--1c05b6e0948ac0c2' }"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar), path(variants)

    output:
    tuple val(meta), path("*.pgen")    , emit: extract_pgen
    tuple val(meta), path("*.psam")    , emit: extract_psam
    tuple val(meta), path("*.pvar.zst"), emit: extract_pvar
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$pgen" == "${prefix}.pgen" ) error "Input and output names are the same, use \"task.ext.prefix\" in modules.config to disambiguate!"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --pfile ${pgen.baseName} \\
        $args \\
        --extract $variants \\
        --make-pgen vzs \\
        --out ${prefix}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst
    """
}
