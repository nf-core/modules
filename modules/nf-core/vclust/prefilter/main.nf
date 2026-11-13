// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VCLUST_PREFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04cb9a3b37b513bfdabb8f0e87050c5573bb7be92926617294cb3de221d37301/data'
        : 'community.wave.seqera.io/library/vclust:1.3.1--039ae9c4872009ec'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('vclust'), eval("vclust --version | sed 's/v//'"), topic: versions, emit: versions_vclust

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vclust \\
        prefilter \\
        $args \\
        -t $task.cpus \\
        -i ${fasta} \\
        -o ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
