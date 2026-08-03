process FGUMI_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64e8f594b6f0dd879bc5abbe4ca70b6b761e1920e407d9e1c7d27b89004aac34/data':
        'community.wave.seqera.io/library/fgumi:0.5.0--a2d14bf52f73eaef' }"

    input:
    tuple val(meta), path(bams, stageAs: "?/*")

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi \\
        merge \\
        --output ${prefix}.bam \\
        --threads $task.cpus \\
        ${args} \\
        ${bams}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
