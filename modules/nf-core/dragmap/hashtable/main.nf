// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DRAGMAP_HASHTABLE {
    tag "${fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // WARN: Do not update this tool to 1.3.0 until https://github.com/Illumina/DRAGMAP/issues/47 is resolved
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e2913ec4630c51abd5e520a2f6c512720d8615f553ad59edc88c90e58df5283/data'
        : 'community.wave.seqera.io/library/dragmap:1.2.1--5d9f1388ce4b8f33'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("dragmap"), emit: hashmap
    tuple val("${task.process}"), val('dragmap'), eval("dragen-os --version 2>&1"), emit: versions_dragmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir dragmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference ${fasta} \\
        --output-directory dragmap \\
        ${args} \\
        --ht-num-threads ${task.cpus}
    """

    stub:
    """
    mkdir dragmap
    """
}
