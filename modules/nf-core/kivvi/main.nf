process KIVVI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d49957352de57e4a080c95d918f16661e1c455def0f7e6c1bc131b32c759e09a/data':
        'community.wave.seqera.io/library/kivvi:1.1.0--d0e9afdd4e10e15c' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val command

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('kivvi'), eval("kivvi --version | sed 's/.* //g'"), topic: versions, emit: versions_kivvi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = command == "d4z4" ? "--threads $task.cpus" : ""
    """
    kivvi \\
        $args \\
        --prefix ${prefix} \\
        --out ./ \\
        --bam $bam \\
        $command \\
        $threads
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
