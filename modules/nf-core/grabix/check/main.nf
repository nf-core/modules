// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GRABIX_CHECK {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f28fcf385db87c5850651ac4600bdd85a749bf67fd55e672541c632ab51d648d/data'
        : 'community.wave.seqera.io/library/grabix:0.1.8--10b827349d23cb56'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), stdout, emit: compress_bgzip
    tuple val("${task.process}"), val('grabix'), eval("grabix |& sed -n 's/^version: //p'"), topic: versions, emit: versions_grabix

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    grabix check ${input} | tr -d '\\n'
    """

    stub:

    """
    echo yes | tr -d '\\n'
    """
}
