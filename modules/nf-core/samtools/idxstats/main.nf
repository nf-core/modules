// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_IDXSTATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        idxstats \\
        --threads ${task.cpus - 1} \\
        ${bam} \\
        > ${prefix}.idxstats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.idxstats
    """
}
