// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLINK_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/5901c2e9a6b11094536c8b0ab0932bd61bb55dbfe9ad270d3babc8f7246a4759/data'
:         'community.wave.seqera.io/library/plink:1.90b6.21--0baa44af6db8f28d' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.genome"), emit: genome
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --genome \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}
    """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   touch ${prefix}.genome
   """

}
