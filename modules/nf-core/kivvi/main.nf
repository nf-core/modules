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
    tuple val(meta), path("*.kivvi.*.json")   , emit: json
    tuple val(meta), path("*.kivvi.*.vcf")    , emit: vcf
    tuple val(meta), path("*.kivvi.*.svg")    , emit: svg, optional: true
    tuple val(meta), path("*.kivvi.*.bam")    , emit: bam
    tuple val(meta), path("*.kivvi.*.bam.bai"), emit: bai
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
    touch ${prefix}.kivvi.${command}.json
    touch ${prefix}.kivvi.${command}.vcf
    touch ${prefix}.kivvi.${command}.svg
    touch ${prefix}.kivvi.${command}.bam
    touch ${prefix}.kivvi.${command}.bam.bai
    """
}
