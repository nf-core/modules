// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTTREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.1.10--h516909a_4' :
        'quay.io/biocontainers/fasttree:2.1.10--h516909a_4' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.nwk"), emit: phylogeny
    tuple val("${task.process}"), val('fasttree'), eval('fasttree -help 2>&1 | head -1 | sed \'s/^FastTree \\([0-9.]*\\) .*$/\\1/\''), topic: versions, emit: versions_fasttree
    tuple val(meta), path("*.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    fasttree \\
        $args \\
        -log ${prefix}_tree.log \\
        -nt $alignment \\
        > ${prefix}.nwk
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"
    """
     touch ${prefix}.nwk
     touch ${prefix}_tree.log
    """
}
