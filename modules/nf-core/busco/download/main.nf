process BUSCO_DOWNLOAD {
    tag "${lineage}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd93359cce52c80ad33d58e18a694f37dc8aee9ccf2530cce8c42f9ade5c1f3f/data'
        : 'community.wave.seqera.io/library/busco:6.1.0--c502c35c2ab05540' }"

    input:
    val lineage

    output:
    path "busco_downloads", emit: download_dir
    tuple val("${task.process}"), val('busco'), eval("busco --version 2> /dev/null | sed 's/BUSCO //g'"), emit: versions_busco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        --download ${lineage} \\
        ${args}
    """

    stub:
    """
    mkdir busco_downloads
    """
}
