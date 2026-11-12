// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KAIJU_MKFMI {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/40cd79a6bb3b5feabe877bef8498d63a3b91fed44b30556af4c7e6554c4310f2/data'
        : 'community.wave.seqera.io/library/kaiju:1.10.3--8b4a2069cc29ac7a'}"

    input:
    tuple val(meta), path(fasta)
    path nodes_dmp, stageAs: "nodes.dmp"
    path names_dmp, stageAs: "names.dmp"
    val keep_intermediate

    output:
    tuple val(meta), path("*.{fmi,dmp}", includeInputs: true), emit: fmi
    tuple val(meta), path("*.bwt"), optional: true, emit: bwt
    tuple val(meta), path("*.sa"), optional: true, emit: sa
    tuple val("${task.process}"), val('kaiju'), eval("kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'"), emit: versions_kaiju, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_cleanup = keep_intermediate ? "" : "rm -f *.{bwt,sa}"
    """
    kaiju-mkbwt \\
        ${args} \\
        -n ${task.cpus} \\
        -o ${prefix} \\
        ${fasta}

    kaiju-mkfmi ${prefix}

    ${run_cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_cleanup = keep_intermediate ? "" : "rm -f *.{bwt,sa}"
    """
    touch ${prefix}.fmi
    touch ${prefix}.bwt
    touch ${prefix}.sa

    ${run_cleanup}
    """
}
