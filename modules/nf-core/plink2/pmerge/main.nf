// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLINK2_PMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dfbbefb7c4dccd8c211a1c8bf94988dcd6a58a147e5c22e2553e80af907569/data'
:         'community.wave.seqera.io/library/plink2:2.00a5.10--223badfc86e2b6e6' }"

    input:
    tuple val(meta), path(pgen, stageAs: 'input/*'), path(pvar, stageAs: 'input/*'), path(psam, stageAs: 'input/*')
    tuple val(meta2), path(pgen2, stageAs: 'merge/*'), path(pvar2, stageAs: 'merge/*'), path(psam2, stageAs: 'merge/*')

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.pvar"), emit: pvar
    tuple val(meta), path("*.psam"), emit: psam
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pgen ${pgen} \\
        --pvar ${pvar} \\
        --psam ${psam} \\
        --pmerge ${pgen2} ${pvar2} ${psam2} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pgen
    touch ${prefix}.pvar
    touch ${prefix}.psam
    """
}
