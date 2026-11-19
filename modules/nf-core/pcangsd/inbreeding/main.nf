// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PCANGSD_INBREEDING {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/853e9d09d2da49363420f6c4a9bc1617145b73ee9dd651add902e5fb5fbdbb9d/data'
        : 'community.wave.seqera.io/library/pcangsd:1.36.4--ef9ce72a07cdfbbd'}"

    input:
    tuple val(meta), path(beagle_file)

    output:
    tuple val(meta), path("*.inbreed.samples"), emit: inbreeding_coefficients
    tuple val("${task.process}"), val('pcangsd'), eval("pcangsd --version"), topic: versions, emit: versions_pcangsd

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pcangsd \\
        --threads ${task.cpus} \\
        --beagle ${beagle_file} \\
        --inbreed-samples \\
        --out ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.inbreed.samples
    """
}
